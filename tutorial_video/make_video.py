"""Build the elastic_vasp tutorial video (fcc Al, real VASP data).

    python3 make_video.py --data /path/to/vasp_alp [--out build] \
        --tts-python /path/to/venv/bin/python --tts-models /path/to/kokoro

The narration uses the Kokoro neural voice (tts_kokoro.py, run with --tts-python, an
interpreter that has kokoro-onnx and soundfile) or, with --engine say, the macOS `say` voices.

The data folder holds the finished calculations:
    std/{C11_C12_I,C11_C12_II,C44}/Result.txt, std/POSCAR, std/post2.out   (cubic_post)
    rot/POSCAR, rot/post2.out                                              (auto_post)
    ev/v*/{POSCAR,OUTCAR}                                                  (energy-volume scan)
    inputs/{INCAR,KPOINTS}

Output (in --out): tutorial.mp4, tutorial.srt, thumbnail.png, youtube_description.txt.
Needs Kokoro (or macOS `say`), ffmpeg, matplotlib, Pillow and ASE.
"""
import argparse
import functools
import io
import json
import os
import re
import subprocess
import sys
import wave
import zlib

import numpy as np
from PIL import Image, ImageDraw, ImageFont

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

W, H, FPS = 1920, 1080, 30
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SF = "/System/Library/Fonts/SFNS.ttf"
MENLO = "/System/Library/Fonts/Menlo.ttc"

BG1, BG2 = (9, 13, 26), (19, 28, 52)
FG = (233, 239, 248)
FG2 = (196, 206, 224)
MUTED = (128, 142, 168)
CYAN = (76, 201, 240)
ORANGE = (255, 159, 67)
GREEN = (87, 204, 153)
RED = (239, 84, 111)
PURPLE = (184, 161, 255)
YELLOW = (255, 209, 102)
PANEL = (21, 29, 51)
BORDER = (50, 63, 96)
TERM_BG = (12, 16, 27)

LEAD, GAP, TAIL = 0.7, 0.45, 1.3
CPS = 30.0  # typing speed, characters per second


def hexc(c):
    return "#%02x%02x%02x" % c


# --------------------------------------------------------------------------- fonts and sprites

@functools.lru_cache(None)
def sf(size, weight="Regular"):
    f = ImageFont.truetype(SF, size)
    for name in f.get_variation_names():
        if (name.decode() if isinstance(name, bytes) else name) == weight:
            f.set_variation_by_name(name)
            break
    return f


@functools.lru_cache(None)
def mono(size, bold=False):
    return ImageFont.truetype(MENLO, size, index=1 if bold else 0)


def ease(x):
    x = min(max(x, 0.0), 1.0)
    return x * x * (3 - 2 * x)


def appear(since, dur=0.45):
    """(opacity, y offset) for an element that starts appearing `since` seconds ago."""
    if since is None or since < 0:
        return 0.0, 0
    p = ease(since / dur)
    return p, int((1 - p) * 22)


_fade_cache = {}


def faded(sprite, op):
    q = int(round(op * 20))
    if q >= 20:
        return sprite
    key = (id(sprite), q)
    if key not in _fade_cache:
        s = sprite.copy()
        a = s.getchannel("A").point([int(v * q / 20) for v in range(256)])
        s.putalpha(a)
        _fade_cache[key] = (sprite, s)
    return _fade_cache[key][1]


def put(cv, sprite, x, y, op=1.0, dy=0):
    if op <= 0.01:
        return
    cv.alpha_composite(faded(sprite, op), (int(x), int(y + dy)))


def put_in(cv, sprite, x, y, since, dur=0.45):
    op, dy = appear(since, dur)
    put(cv, sprite, x, y, op, dy)


def wrap(text, f, maxw):
    if maxw is None:
        return text.split("\n")
    out = []
    for para in text.split("\n"):
        line = ""
        for word in para.split(" "):
            trial = (line + " " + word).strip()
            if f.getlength(trial) <= maxw or not line:
                line = trial
            else:
                out.append(line)
                line = word
        out.append(line)
    return out


@functools.lru_cache(None)
def text_sprite(text, size, color=FG, weight="Regular", is_mono=False, maxw=None, spacing=1.32, align="left"):
    f = mono(size, weight == "Bold") if is_mono else sf(size, weight)
    lines = wrap(text, f, maxw)
    lh = int(size * spacing)
    w = int(max(f.getlength(line) for line in lines)) + 4
    if maxw and align == "center":
        w = maxw
    h = lh * len(lines) + int(size * 0.35)
    img = Image.new("RGBA", (max(w, 1), h), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    for i, line in enumerate(lines):
        x = (w - f.getlength(line)) / 2 if align == "center" else 0
        d.text((x, i * lh), line, font=f, fill=color)
    return img


def rich_sprite(parts, size, gap=0):
    """One line made of (text, color, weight[, mono]) parts."""
    imgs = [text_sprite(p[0], size, p[1], p[2] if len(p) > 2 else "Regular", len(p) > 3 and p[3]) for p in parts]
    w = sum(i.width for i in imgs) + gap * (len(imgs) - 1)
    h = max(i.height for i in imgs)
    img = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    x = 0
    for i in imgs:
        img.alpha_composite(i, (x, h - i.height))
        x += i.width + gap
    return img


@functools.lru_cache(None)
def chip(text, size=28, fg=FG, bg=PANEL, border=BORDER, weight="Medium", is_mono=False, pad=(22, 10)):
    t = text_sprite(text, size, fg, weight, is_mono, spacing=1.15)
    w, h = t.width + 2 * pad[0], t.height + 2 * pad[1] - int(size * 0.3)
    img = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    d.rounded_rectangle((0, 0, w - 1, h - 1), radius=h // 2, fill=bg + (255,), outline=border, width=2)
    img.alpha_composite(t, (pad[0], pad[1] - 2))
    return img


def panel_img(w, h, fill=PANEL, border=BORDER, radius=20, width=2):
    img = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    ImageDraw.Draw(img).rounded_rectangle((0, 0, w - 1, h - 1), radius=radius, fill=fill + (255,),
                                          outline=border, width=width)
    return img


@functools.lru_cache(None)
def code_panel(title, text, w, size=27, hl=(), hl_color=YELLOW, h=None):
    lines = text.split("\n")
    lh = int(size * 1.48)
    hh = h or 66 + lh * len(lines) + 24
    img = panel_img(w, hh, fill=TERM_BG)
    d = ImageDraw.Draw(img, "RGBA")
    d.text((26, 16), title, font=sf(24, "Semibold"), fill=MUTED)
    d.line((0, 56, w, 56), fill=BORDER, width=2)
    for i, line in enumerate(lines):
        y = 70 + i * lh
        if i in hl:
            d.rounded_rectangle((12, y - 6, w - 12, y + lh - 6), radius=8, fill=hl_color + (40,))
            d.rectangle((12, y - 6, 17, y + lh - 6), fill=hl_color + (255,))
        comment = line.find("!")
        if comment >= 0:
            d.text((30, y), line[:comment], font=mono(size), fill=FG)
            d.text((30 + mono(size).getlength(line[:comment]), y), line[comment:], font=mono(size), fill=MUTED)
        else:
            d.text((30, y), line, font=mono(size), fill=FG)
    return img


@functools.lru_cache(None)
def math_sprite(tex, size, color=FG):
    plt.rcParams["mathtext.fontset"] = "cm"
    fig = plt.figure(figsize=(0.01, 0.01), dpi=100)
    fig.text(0, 0, tex, fontsize=size * 0.72, color=hexc(color))
    buf = io.BytesIO()
    fig.savefig(buf, dpi=100, transparent=True, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGBA")


@functools.lru_cache(None)
def matrix_sprite(rows, size=30, colors=None, colw=None, bracket=FG2):
    f = sf(size, "Regular")
    ncol = len(rows[0])
    widths = [colw or max(f.getlength(r[j]) for r in rows) + size * 1.0 for j in range(ncol)]
    rh = int(size * 1.55)
    w = int(sum(widths)) + 40
    h = rh * len(rows) + 16
    img = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    for i, r in enumerate(rows):
        x = 20
        for j, cell in enumerate(r):
            col = colors[i][j] if colors else FG
            d.text((x + (widths[j] - f.getlength(cell)) / 2, 8 + i * rh + rh * 0.12), cell, font=f, fill=col)
            x += widths[j]
    for x0, sgn in ((4, 1), (w - 5, -1)):
        d.line([(x0 + sgn * 12, 4), (x0, 4), (x0, h - 5), (x0 + sgn * 12, h - 5)], fill=bracket, width=3)
    return img


def image_card(path, maxw, maxh, pad=14):
    im = Image.open(path).convert("RGBA")
    s = min((maxw - 2 * pad) / im.width, (maxh - 2 * pad) / im.height)
    im = im.resize((int(im.width * s), int(im.height * s)), Image.LANCZOS)
    card = panel_img(im.width + 2 * pad, im.height + 2 * pad, fill=(255, 255, 255), border=(255, 255, 255), radius=18)
    card.alpha_composite(im, (pad, pad))
    return card


# --------------------------------------------------------------------------- terminal

class Term:
    def __init__(self, x, y, w, h, title="bash — Al", size=25):
        self.x, self.y, self.w, self.h, self.title, self.size = x, y, w, h, title, size
        self.lh = int(size * 1.5)
        self.items = []
        self.cache = {}

    def cmd(self, text, t, prompt="$ "):
        self.items.append(dict(k="cmd", s=text, t=t, p=prompt))
        return t + len(text) / CPS

    def out(self, lines, t, stagger=0.05):
        for i, line in enumerate(lines):
            if isinstance(line, str):
                line = dict(s=line)
            self.items.append(dict(k="out", s=line["s"], t=t + i * stagger, c=line.get("c", FG2),
                                   hl=line.get("hl")))
        return t + len(lines) * stagger

    def state(self, t):
        rows, typing = [], False
        for it in self.items:
            if t < it["t"]:
                continue
            if it["k"] == "cmd":
                n = int((t - it["t"]) * CPS)
                rows.append((it["p"], it["s"][:n], FG, None))
                if n < len(it["s"]):
                    typing = True
                    break
            else:
                hl = it["hl"][1] if it["hl"] and t >= it["hl"][0] else None
                rows.append(("", it["s"], it["c"], hl))
        cursor = typing or (int(t * 2) % 2 == 0)
        if not typing:
            rows.append(("$ ", "", FG, None))
        maxrows = (self.h - 70) // self.lh
        return tuple(rows[-maxrows:]), cursor

    def sprite(self, t):
        key = self.state(t)
        if key not in self.cache:
            self.cache[key] = self._render(*key)
        return self.cache[key]

    def _render(self, rows, cursor):
        img = panel_img(self.w, self.h, fill=TERM_BG, border=BORDER, radius=16)
        d = ImageDraw.Draw(img, "RGBA")
        for i, c in enumerate(((255, 95, 87), (254, 188, 46), (40, 200, 64))):
            d.ellipse((20 + i * 26, 16, 34 + i * 26, 30), fill=c)
        tw = sf(20).getlength(self.title)
        d.text(((self.w - tw) / 2, 12), self.title, font=sf(20), fill=MUTED)
        d.line((0, 46, self.w, 46), fill=BORDER, width=2)
        f = mono(self.size)
        cw = f.getlength("M")
        for i, (prefix, text, color, hl) in enumerate(rows):
            y = 60 + i * self.lh
            if hl:
                d.rounded_rectangle((10, y - 5, self.w - 10, y + self.lh - 7), radius=6, fill=hl + (46,))
                d.rectangle((10, y - 5, 14, y + self.lh - 7), fill=hl + (255,))
            x = 24
            if prefix:
                d.text((x, y), prefix, font=mono(self.size, True), fill=GREEN)
                x += cw * len(prefix)
            d.text((x, y), text, font=f, fill=color)
            if i == len(rows) - 1 and prefix and cursor:
                cx = x + cw * len(text)
                d.rectangle((cx + 2, y + 2, cx + cw, y + self.size + 4), fill=FG + (220,))
        return img

    def draw(self, cv, t, since=0.0):
        op, dy = appear(since)
        put(cv, self.sprite(t), self.x, self.y, op, dy)


# --------------------------------------------------------------------------- plots

def style_ax(ax):
    ax.set_facecolor("none")
    for s in ax.spines.values():
        s.set_color(hexc(BORDER))
    ax.tick_params(colors=hexc(FG2), labelsize=13)
    ax.xaxis.label.set_color(hexc(FG2))
    ax.yaxis.label.set_color(hexc(FG2))
    ax.grid(True, color=hexc(MUTED), alpha=0.18, lw=0.8)


def fig_to_sprite(fig):
    buf = io.BytesIO()
    fig.savefig(buf, dpi=100, transparent=True)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGBA")


def strain_plot(ax, d, e, npts, frac, color, label=None):
    style_ax(ax)
    ax.set_xlim(-3.4, 3.4)
    ax.set_ylim(-0.08 * e.max(), e.max() * 1.12)
    p = np.polyfit(d, e, 2)
    if frac > 0:
        x = np.linspace(-3.3, -3.3 + 6.6 * frac, 200)
        ax.plot(x, np.polyval(p, x), color=hexc(color), lw=2.6, alpha=0.95)
    if npts:
        order = np.argsort(np.abs(d))  # from the centre outwards
        idx = order[:npts]
        ax.plot(d[idx], e[idx], "o", ms=10, color=hexc(FG), mec=hexc(color), mew=2.5, zorder=5)
    ax.set_xlabel("strain δ (%)", fontsize=15)
    ax.set_ylabel("E − E₀ (meV)", fontsize=15)
    if label:
        ax.set_title(label, color=hexc(FG), fontsize=17, loc="left", pad=10)


# --------------------------------------------------------------------------- data

def read_data(root):
    from ase.io import read
    from ase.eos import EquationOfState

    data = {}
    curves = {}
    for name in ("C11_C12_I", "C11_C12_II", "C44"):
        arr = np.loadtxt(os.path.join(root, "std", name, "Result.txt"))
        curves[name] = (arr[:, 0] * 100, arr[:, 1] * 1000)  # %, meV
    data["curves"] = curves

    def numbers(path):
        txt = open(path).read()
        return {k + sl: float(v) for k, sl, v in re.findall(r"^(\S+?)(/?)= *([-\d.]+)", txt, re.M)}, txt

    data["std"], data["std_txt"] = numbers(os.path.join(root, "std", "post2.out"))
    rot, rot_txt = numbers(os.path.join(root, "rot", "post2.out"))
    data["rot"] = rot
    data["rot_txt"] = rot_txt

    vols, ens = [], []
    for v in sorted(os.listdir(os.path.join(root, "ev"))):
        folder = os.path.join(root, "ev", v)
        vols.append(read(os.path.join(folder, "POSCAR"), format="vasp").get_volume())
        ens.append(float(re.findall(r"energy\(sigma->0\) =\s*([-\d.]+)", open(os.path.join(folder, "OUTCAR")).read())[-1]))
    eos = EquationOfState(vols, ens, eos="birchmurnaghan")
    v0, e0, b = eos.fit()
    data["eos"] = dict(v=np.array(vols), e=np.array(ens), v0=v0, e0=e0, B=b * 160.21766, a0=(4 * v0) ** (1 / 3), fit=eos)
    lines = open(os.path.join(root, "std", "POSCAR")).read().rstrip("\n").split("\n")
    for i in (2, 3, 4, 8):
        lines[i] = "  " + "  ".join("%10.6f" % (float(v) + 0.0) for v in lines[i].split()[:3])
    lines[1] = "  " + "%.6f" % float(lines[1])
    data["poscar_std"] = "\n".join(line.rstrip() for line in lines)
    data["lat_std"] = read(os.path.join(root, "std", "POSCAR"), format="vasp").cell.array
    data["lat_rot"] = read(os.path.join(root, "rot", "POSCAR"), format="vasp").cell.array
    short = {"ISMEAR": "! tetrahedron method", "NSW": "! 1 atom: no relaxation"}
    incar = []
    for line in open(os.path.join(root, "inputs", "INCAR")).read().splitlines():
        if not line.strip() or line.startswith(("SYSTEM", "LSCALAPACK", "LWAVE", "LCHARG")):
            continue
        line = line.split("!")[0].rstrip()
        key = line.split("=")[0].strip()
        incar.append(line.ljust(16) + short[key] if key in short else line)
    data["incar"] = "\n".join(incar)
    data["kpoints"] = open(os.path.join(root, "inputs", "KPOINTS")).read().rstrip("\n")
    data["check"] = open(os.path.join(root, "check_std.out")).read().replace("std/POSCAR", "POSCAR").rstrip("\n").split("\n")
    return data


# --------------------------------------------------------------------------- speech

DIGITS = dict(zip("123456", ("one", "two", "three", "four", "five", "six")))
SPOKEN = [
    (r"y = aδ² \+ bδ \+ c", "a delta squared, plus b delta, plus c"),
    (r"elastic_vasp", "elastic vasp"), (r"elastic-vasp", "elastic vasp"), (r"elastic-vis", "elastic vis"),
    (r"\bGPa\b", "gigapascals"), (r"Fm-3m", "F m bar 3 m"), (r"spglib", "S P G lib"), (r"\bASE\b", "A S E"),
    (r"INCAR", "in-car"), (r"KPOINTS", "k-points"), (r"POTCAR", "pot-car"), (r"OUTCAR", "out-car"),
    (r"POSCAR", "pos-car"), (r"\bVASP\b", "vasp"), (r"δ", "delta"), (r"\bC([1-6])([1-6])\b", lambda m: "C %s %s" % (DIGITS[m[1]], DIGITS[m[2]])),
    (r"a₀", "a zero"), (r"Å", "angstrom"), (r"\beV\b", "electron-volt"),
    (r"(\d+)×\1×\1", r"\1 by \1 by \1"), (r"ISMEAR = −5", "ISMEAR minus five"), (r"ISIF = 2", "ISIF two"),
    (r"Birch–Murnaghan", "Birch Murnegan"), (r"\borthorhombic\b", "ortho-rombic"),
    (r"\btrigonal\b", "trygonal"),
    (r"\beigenvalue", "eigan value"), (r"Γ", "gamma"), (r"--output", "the output option"), (r"--help", "with the help option"), (r"–", " "), (r"Python 3\.10", "Python three ten"),
    (r"(\d+)\.0\b", r"\1"), (r" = ", " equals "),
]


def spoken(text):
    for a, b in SPOKEN:
        text = re.sub(a, b, text)
    return text


def synth_all(jobs, args):
    """Write the missing clips of ``jobs`` ([(text, wav_path)]) with the chosen engine."""
    jobs = [j for j in jobs if not os.path.exists(j[1])]
    if not jobs:
        return
    if args.engine == "say":
        for text, path in jobs:
            subprocess.run(["say", "-v", args.voice, "-r", str(args.rate), "-o", path, "--file-format=WAVE",
                            "--data-format=LEI16@24000", text], check=True)
        return
    job_file = os.path.join(args.out, "voice", "jobs.json")
    with open(job_file, "w") as f:
        json.dump(jobs, f)
    here = os.path.dirname(os.path.abspath(__file__))
    if args.engine == "chatterbox":
        subprocess.run([args.tts_python, os.path.join(here, "tts_chatterbox.py"), args.voice_sample,
                        str(args.exaggeration), str(args.cfg), job_file], check=True)
    else:
        subprocess.run([args.tts_python, os.path.join(here, "tts_kokoro.py"), args.tts_models, args.voice,
                        str(args.speed), job_file], check=True)


SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "narration.txt")


def apply_script(scenes, path=SCRIPT):
    """Use the sentences of the editable script ``path`` (written from the code if missing).

    A line that differs from the code replaces that sentence (subtitle and narration); the
    pronunciation rules in SPOKEN still apply to it."""
    if not os.path.exists(path):
        write_script(path, scenes)
        return
    text = {}
    for line in open(path, encoding="utf-8"):
        m = re.match(r"(s\d\d_\d+)\s+(.*\S)", line)
        if m:
            text[m[1]] = re.sub(r"\s+", " ", m[2])
    for k, sc in enumerate(scenes):
        sc.sentences = [(text[key], None, pause) if text.get(key, sub) != sub else (sub, say, pause)
                        for j, (sub, say, pause) in enumerate(sc.sentences)
                        for key in ["s%02d_%d" % (k + 1, j + 1)]]


def write_script(path, scenes):
    """Numbered narration script: edit the text after an id to change that sentence."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("Narration of the tutorial video. Edit the text after an id to change that sentence,\n"
                "then run make_video.py again (this file is never overwritten).\n"
                "To use your own voice, record each line as a file named by its id (e.g. s02_1.m4a)\n"
                "into one folder and run make_video.py with --recorded FOLDER.\n")
        for k, sc in enumerate(scenes):
            f.write("\n== %d. %s ==\n" % (k + 1, sc.chapter))
            for j, (sub, say_text, pause) in enumerate(sc.sentences):
                f.write("s%02d_%d   %s\n" % (k + 1, j + 1, sub))


def recorded_clip(args, k, j):
    """Your own recording of sentence j of scene k, converted to 24 kHz mono WAVE (None if absent)."""
    for ext in ("wav", "m4a", "mp3", "aiff", "aif", "caf", "flac"):
        src = os.path.join(args.recorded, "s%02d_%d.%s" % (k + 1, j + 1, ext))
        if os.path.exists(src):
            dst = os.path.join(args.out, "voice", "rec_s%02d_%d.wav" % (k + 1, j + 1))
            if not os.path.exists(dst) or os.path.getmtime(dst) < os.path.getmtime(src):
                trim = "silenceremove=start_periods=1:start_threshold=-45dB,areverse," \
                       "silenceremove=start_periods=1:start_threshold=-45dB,areverse"
                tmp = dst[:-4] + "_raw.wav"   # macOS afconvert reads every Apple/QuickTime format
                subprocess.run(["afconvert", "-f", "WAVE", "-d", "LEI16@24000", "-c", "1", src, tmp], check=True)
                subprocess.run([args.ffmpeg, "-v", "error", "-y", "-i", tmp, "-af", trim, "-ac", "1", "-ar", "24000",
                                "-sample_fmt", "s16", dst], check=True)
                os.remove(tmp)
            return dst
    return None


def read_clip(path):
    with wave.open(path) as w:
        assert w.getframerate() == 24000 and w.getsampwidth() == 2 and w.getnchannels() == 1, path
        return w.readframes(w.getnframes()), w.getnframes() / w.getframerate()


# --------------------------------------------------------------------------- scenes

class Scene:
    chapter = ""
    title = ""
    sentences = []   # (subtitle, spoken override or None, extra pause after)

    def setup(self, D):
        pass

    def c(self, i, plus=0.0):
        return self.t - self.cues[i] - plus

    def draw(self, cv, D):
        pass


class Intro(Scene):
    chapter = "Introduction"
    title = ""
    sentences = [
        ("In this tutorial we compute the elastic constants of a crystal with VASP and elastic_vasp.", None, 0.3),
        ("Every number you will see comes from real VASP calculations on fcc aluminium.", None, 0.6),
    ]

    def setup(self, D):
        self.t1 = text_sprite("Elastic constants\nwith VASP", 104, FG, "Bold", spacing=1.12)
        self.t2 = rich_sprite([("elastic_vasp", CYAN, "Bold", True), ("  ·  a step-by-step tutorial", FG2, "Medium")], 40)
        self.t3 = text_sprite("Worked example: fcc aluminium — real VASP results", 32, MUTED, "Regular")
        self.img = image_card(os.path.join(REPO, "docs/_static/elastic_vis/Al_Youngs_Modulus_3D.png"), 720, 720)
        self.chips = [chip("check", 30, CYAN, is_mono=True), chip("pre", 30, ORANGE, is_mono=True),
                      chip("VASP", 30, FG, is_mono=True), chip("post", 30, GREEN, is_mono=True)]

    def draw(self, cv, D):
        put_in(cv, self.t1, 130, 260, self.t - 0.1, 0.7)
        put_in(cv, self.t2, 134, 545, self.t - 0.5, 0.7)
        put_in(cv, self.t3, 134, 620, self.t - 0.8, 0.7)
        x = 134
        for i, ch in enumerate(self.chips):
            put_in(cv, ch, x, 730, self.c(0, 1.5 + 0.35 * i))
            x += ch.width + 46
            if i < 3:
                put_in(cv, text_sprite("→", 34, MUTED), x - 38, 733, self.c(0, 1.6 + 0.35 * i))
        put_in(cv, self.img, 1110, 200, self.c(1), 0.8)


class Idea(Scene):
    chapter = "The method"
    title = "The idea: energy versus strain"
    sentences = [
        ("Distort the relaxed crystal by a series of small strains δ, and compute the change in energy for each one.", None, 0.2),
        ("The elastic constants are the second derivatives of this strain energy with respect to strain.", None, 0.3),
        ("So for each distortion we fit the energy change per volume with a quadratic, y = aδ² + bδ + c, and keep the curvature a.",
         "So for each distortion we fit the energy change per volume with a quadratic, a delta squared, plus b delta, plus c, and keep the curvature a.", 0.3),
        ("For a cubic crystal we use three independent distortions: orthorhombic, hydrostatic and trigonal shear.", None, 0.3),
        ("Each curvature is a known combination of C11, C12 and C44: solve the three equations and you have the full stiffness matrix.", None, 0.8),
    ]

    def setup(self, D):
        self.m1 = math_sprite(r"$C_{ij} \;=\; \dfrac{1}{V_0}\,\dfrac{\partial^2 E}{\partial \varepsilon_i\,\partial \varepsilon_j}$", 60)
        self.m2 = math_sprite(r"$y(\delta) = \dfrac{E(\delta) - E_0}{V_0} = a\,\delta^2 + b\,\delta + c$", 50, CYAN)
        modes = (("orthorhombic", r"$\varepsilon = (\delta, -\delta, 0, 0, 0, 0)$", r"$a = C_{11} - C_{12}$", CYAN),
                 ("hydrostatic", r"$\varepsilon = (\delta, \delta, \delta, 0, 0, 0)$", r"$a = \frac{3}{2}(C_{11} + 2C_{12})$", ORANGE),
                 ("trigonal shear", r"$\varepsilon = (0, 0, 0, \delta, \delta, \delta)$", r"$a = \frac{3}{2}C_{44}$", PURPLE))
        self.modes = [(text_sprite(lab, 28, col, "Semibold"), math_sprite(eps, 34, FG2), math_sprite(av, 34, FG))
                      for lab, eps, av, col in modes]
        self.d, self.e = D["curves"]["C44"]
        self.cache = {}
        steps = ["1. Strain by δ", "2. Compute E", "3. Fit a", "4. Solve for C"]
        self.steps = [chip(s, 30, FG, PANEL, BORDER, "Semibold") for s in steps]

    def plot(self, n, frac):
        key = (n, round(frac * 14))
        if key not in self.cache:
            fig = plt.figure(figsize=(8.4, 5.6), dpi=100)
            ax = fig.add_axes([0.13, 0.15, 0.82, 0.8])
            strain_plot(ax, self.d, self.e, n, key[1] / 14, ORANGE)
            self.cache[key] = fig_to_sprite(fig)
        return self.cache[key]

    def draw(self, cv, D):
        n = int(np.clip(self.c(0, 0.8) / 0.25 + 1, 0, 7)) if self.c(0, 0.8) >= 0 else 0
        frac = ease(self.c(2, 1.0) / 1.2) if self.c(2, 1.0) > 0 else 0
        put_in(cv, self.plot(n, frac), 1010, 250, self.c(0))
        put_in(cv, self.m1, 140, 235, self.c(1))
        put_in(cv, self.m2, 140, 410, self.c(2))
        for i, (lab, eps, av) in enumerate(self.modes):
            y, since = 560 + i * 92, self.c(3, 2.4 + 0.9 * i)
            put_in(cv, lab, 140, y + 10, since)
            put_in(cv, eps, 360, y + 24 - eps.height // 2, since)
            put_in(cv, av, 720, y + 24 - av.height // 2, since)
        x = 140
        for i, s in enumerate(self.steps):
            since = self.c(0, 0.2) if i == 0 else self.c(0, 1.6) if i == 1 else self.c(2, 1.0) if i == 2 else self.c(4, 2.0)
            put_in(cv, s, x, 885, since)
            x += s.width + 34


class Install(Scene):
    chapter = "Install"
    title = "Install"
    sentences = [
        ("Install the package from GitHub with pip.", None, 1.4),
        ("It needs Python 3.10 or newer, and installs ASE and spglib with it.", None, 0.4),
        ("Then check that the elastic-vasp command is available.", None, 2.2),
    ]

    def setup(self, D):
        T = self.term = Term(110, 240, 1180, 450)
        t = T.cmd('pip install "git+https://github.com/prnvrvs/elastic_vasp.git"', self.cues[0] + 0.3)
        T.out([dict(s="Successfully installed elastic-vasp-0.1.1", c=GREEN)], t + 0.5)
        t = T.cmd("elastic-vasp --help", self.cues[2] + 0.4)
        T.out(["usage: elastic-vasp [-h] [--version] [-C DIRECTORY] COMMAND ...", "",
               "  check    report the crystal symmetry and which workflow to use",
               "  pre      write the strained POSCARs (pre-processing)",
               "  post     fit the VASP energies and print the elastic constants",
               "  systems  list the crystal systems and their deformations"], t + 0.3)
        self.chips = [chip(s, 30, FG, PANEL, BORDER) for s in ("Python ≥ 3.10", "ASE ≥ 3.23", "spglib ≥ 2.4")]
        self.lbl = text_sprite("Requirements", 28, MUTED, "Semibold")

    def draw(self, cv, D):
        self.term.draw(cv, self.t, self.t - 0.1)
        put_in(cv, self.lbl, 1370, 300, self.c(1))
        for i, ch in enumerate(self.chips):
            put_in(cv, ch, 1370, 360 + i * 92, self.c(1, 0.3 + 0.5 * i))


class Structure(Scene):
    chapter = "Step 0 · Structure"
    title = "Step 0 · Start from a relaxed structure"
    sentences = [
        ("Start from a relaxed structure: here, the one-atom primitive cell of aluminium.", None, 0.3),
        ("An energy–volume scan, fitted with a Birch–Murnaghan equation of state, gives a₀ = 4.039 Å and a bulk modulus of 77.0 GPa.", None, 0.3),
        ("Keep that number in mind: we will compare it with the elastic constants at the end.", None, 0.8),
    ]

    def setup(self, D):
        self.poscar = code_panel("POSCAR", D["poscar_std"], 820, 25)
        E = D["eos"]
        self.E = E
        self.cache = {}
        self.badge = chip("B₀ = %.1f GPa  ·  compare later" % E["B"], 30, YELLOW, (45, 40, 20), YELLOW, "Semibold")
        self.cap = text_sprite("a₀ = %.4f Å     V₀ = %.3f Å³" % (E["a0"], E["v0"]), 32, FG, "Medium")

    def plot(self, n, frac):
        key = (n, round(frac * 14))
        if key not in self.cache:
            E = self.E
            fig = plt.figure(figsize=(8.6, 5.4), dpi=100)
            ax = fig.add_axes([0.15, 0.16, 0.8, 0.78])
            style_ax(ax)
            vv = np.linspace(E["v"].min() - 0.15, E["v"].max() + 0.15, 200)
            ax.set_xlim(vv[0], vv[-1])
            ax.set_ylim((E["e"].min() - E["e0"]) * 1000 - 0.6, (E["e"].max() - E["e0"]) * 1000 * 1.15)
            if key[1]:
                m = int(200 * key[1] / 14)
                ax.plot(vv[:m], (E["fit"].func(vv[:m], *E["fit"].eos_parameters) - E["e0"]) * 1000, color=hexc(CYAN), lw=2.6)
            if n:
                ax.plot(E["v"][:n], (E["e"][:n] - E["e0"]) * 1000, "o", ms=10, color=hexc(FG), mec=hexc(CYAN), mew=2.5)
            ax.set_xlabel("volume (Å³)", fontsize=15)
            ax.set_ylabel("E − E₀ (meV)", fontsize=15)
            self.cache[key] = fig_to_sprite(fig)
        return self.cache[key]

    def draw(self, cv, D):
        put_in(cv, self.poscar, 110, 250, self.c(0))
        n = int(np.clip(self.c(1, 0.2) / 0.2 + 1, 0, 7)) if self.c(1, 0.2) >= 0 else 0
        frac = ease(self.c(1, 1.8) / 1.4) if self.c(1, 1.8) > 0 else 0
        put_in(cv, self.plot(n, frac), 990, 240, self.c(1))
        put_in(cv, self.cap, 1100, 790, self.c(1, 3.5))
        put_in(cv, self.badge, 1100, 870, self.c(2, 0.2))


class Check(Scene):
    chapter = "Step 1 · Check"
    title = "Step 1 · Check the symmetry"
    sentences = [
        ("First, run elastic-vasp check in the folder with your POSCAR.", None, 1.6),
        ("spglib finds the space group, Fm-3m, and the number of independent constants: three for a cubic crystal.", None, 0.3),
        ("It also recommends a workflow, cubic, and confirms that it works in this orientation.", None, 1.0),
    ]

    def setup(self, D):
        T = self.term = Term(110, 240, 1700, 430, size=24)
        t = T.cmd("elastic-vasp check", self.cues[0] + 0.5)
        lines = D["check"]
        out = []
        for i, line in enumerate(lines):
            hl = None
            if i in (1, 2):
                hl = (self.cues[1] + 0.3 + 1.3 * (i - 1), CYAN)
            if i == 3:
                hl = (self.cues[2] + 0.3, GREEN)
            out.append(dict(s=line, c=FG if i else FG2, hl=hl))
        T.out(out, t + 0.4, 0.12)
        self.chips = [chip("Fm-3m (225)", 32, CYAN, is_mono=True), chip("3 constants: C11, C12, C44", 32, FG),
                      chip("→  cubic_pre / cubic_post", 32, GREEN, is_mono=True)]

    def draw(self, cv, D):
        self.term.draw(cv, self.t, self.t - 0.1)
        put_in(cv, self.chips[0], 110, 740, self.c(1, 0.6))
        put_in(cv, self.chips[1], 110 + self.chips[0].width + 30, 740, self.c(1, 2.2))
        put_in(cv, self.chips[2], 110, 840, self.c(2, 0.6))


class Generate(Scene):
    chapter = "Step 2 · Generate"
    title = "Step 2 · Generate the strained structures"
    sentences = [
        ("Now generate the strained structures with elastic-vasp pre cubic.", None, 1.2),
        ("We use strains up to three percent, at seven points: three deformations times seven strains, so twenty-one VASP runs.", None, 0.3),
        ("The strain ε is the small-strain tensor in Voigt notation, with engineering shear strains: ε₄ = 2ε_yz, and so on.",
         "The strain epsilon is the small strain tensor, written in Voigt notation, with engineering shear strains: epsilon four is two times epsilon y z, and so on.", 0.3),
        ("The deformation matrix is F = I + ε, so its off-diagonal entries are half the shear strains, and each POSCAR gets the lattice vectors L times F.",
         "The deformation matrix F is the identity plus the strain tensor, so its off-diagonal entries are half the shear strains, and each pos-car gets the lattice vectors L times F.", 0.3),
        ("One deformation changes the volume; the other two change the shape.", None, 1.0),
    ]

    def setup(self, D):
        T = self.term = Term(110, 225, 1700, 190, size=25)
        t = T.cmd("elastic-vasp pre cubic --max-strain 0.03 --points 7", self.cues[0] + 0.4)
        T.out([dict(s="Deformation folders written. Run VASP in every strain_* folder, then: elastic-vasp post cubic", c=GREEN)], t + 0.4)
        tree = ("Al/\n├── POSCAR\n├── C11_C12_I/\n│   ├── strain_-0.03/   POSCAR\n│   ├── strain_-0.02/\n"
                "│   ├──   …\n│   └── strain_0.03/\n├── C11_C12_II/   7 strains\n└── C44/          7 strains")
        self.tree = code_panel("21 VASP runs", tree, 600, 25)
        d = "δ"
        mats = [
            ("C11_C12_I", (("1+" + d, "0", "0"), ("0", "1−" + d, "0"), ("0", "0", "1")), "shape", CYAN),
            ("C11_C12_II", (("1+" + d, "0", "0"), ("0", "1+" + d, "0"), ("0", "0", "1+" + d)), "volume", ORANGE),
            ("C44", (("1", d + "/2", d + "/2"), (d + "/2", "1", d + "/2"), (d + "/2", d + "/2", "1")), "shape", CYAN),
        ]
        self.mats = []
        for name, rows, kind, col in mats:
            colors = tuple(tuple(MUTED if c == "0" else FG for c in r) for r in rows)
            self.mats.append((text_sprite(name, 26, FG2, "Bold", True), matrix_sprite(rows, 30, colors, colw=92),
                              chip(kind, 26, col, PANEL, col, "Semibold")))
        self.F = math_sprite(r"$F = I + \varepsilon, \qquad L' = L\,F$", 44, FG)
        self.voigt = math_sprite(r"$\varepsilon = (\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, 2\varepsilon_{yz}, "
                                 r"2\varepsilon_{xz}, 2\varepsilon_{xy})$", 32, FG2)

    def draw(self, cv, D):
        self.term.draw(cv, self.t, self.t - 0.1)
        put_in(cv, self.tree, 110, 455, self.c(1, 0.5))
        put_in(cv, self.voigt, 780, 900, self.c(2, 0.5))
        put_in(cv, self.F, 780, 465, self.c(3, 0.2))
        for i, (lbl, m, tag) in enumerate(self.mats):
            x = 780 + i * 355
            put_in(cv, lbl, x, 570, self.c(3, 1.0 + 0.6 * i))
            put_in(cv, m, x - 10, 620, self.c(3, 1.0 + 0.6 * i))
            put_in(cv, tag, x + 60, 820, self.c(4, 0.3 if i == 1 else 2.0))


class Deform(Scene):
    chapter = "How F changes the cell"
    title = "How F changes the shape of the cell"
    sentences = [
        ("Here is what each deformation matrix does to the cell, exaggerated so that you can see it.", None, 0.4),
        ("The orthorhombic distortion stretches the cell along x and shrinks it along y: the volume stays nearly the same, only the shape changes.", None, 0.4),
        ("The hydrostatic distortion scales all three axes equally: the cell stays cubic, and only its volume changes.", None, 0.4),
        ("The trigonal shear tilts the cube along its body diagonal: the right angles between the edges open up or close, and the cube becomes a rhombohedron.", None, 0.4),
        ("The atoms keep their fractional coordinates, so they move with the cell.", None, 1.0),
    ]
    MODES = (("orthorhombic", CYAN, 0.25, lambda d: np.diag([1 + d, 1 - d, 1.0]),
              r"$F = \mathrm{diag}(1{+}\delta,\ 1{-}\delta,\ 1)$", "shape changes, volume ≈ same"),
             ("hydrostatic", ORANGE, 0.18, lambda d: (1 + d) * np.eye(3),
              r"$F = (1{+}\delta)\,I$", "volume changes, shape same"),
             ("trigonal shear", PURPLE, 0.5, lambda d: np.eye(3) + d / 2 * (np.ones((3, 3)) - np.eye(3)),
              r"$F_{ij} = \delta/2\ \ (i \neq j)$", "angles change, cube → rhombohedron"))

    def setup(self, D):
        self.cache = {}
        self.labels = [(text_sprite(name, 32, col, "Semibold"), math_sprite(tex, 34, FG),
                        chip(eff, 26, col, PANEL, col, "Semibold")) for name, col, _, _, tex, eff in self.MODES]
        self.note = text_sprite("Deformations exaggerated (δ up to 0.5); the real runs use |δ| ≤ 0.03.", 26, MUTED)

    def cube(self, i, q):
        key = (i, q)
        if key not in self.cache:
            name, col, dmax, Fd, _, _ = self.MODES[i]
            fig = plt.figure(figsize=(5.4, 5.0), dpi=100)
            ax = fig.add_axes([0, 0, 1, 1], projection="3d")
            ax.set_axis_off()
            ax.view_init(elev=22, azim=-58)
            ax.set_box_aspect((1, 1, 1))
            for lim in (ax.set_xlim, ax.set_ylim, ax.set_zlim):
                lim(-0.78, 0.78)
            corners = np.array([(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)], float)
            edges = [(m, n) for m in range(8) for n in range(m + 1, 8) if np.abs(corners[m] - corners[n]).sum() == 1]
            for F, style in ((np.eye(3), dict(color=hexc(MUTED), lw=1.6, ls="--")),
                             (Fd(dmax * q / 12), dict(color=hexc(col), lw=3.2))):
                v = corners @ F
                v -= v.mean(axis=0)
                for m, n in edges:
                    ax.plot(*zip(v[m], v[n]), **style)
                if style.get("ls") is None:
                    ax.scatter(*v.T, s=90, color=hexc(FG), edgecolors=hexc(col), linewidths=2, depthshade=False)
            self.cache[key] = fig_to_sprite(fig)
        return self.cache[key]

    def draw(self, cv, D):
        for i in range(3):
            since = self.c(i + 1, 0.6)
            q = int(round(12 * ease(min(1.0, since / 1.8)))) if since > 0 else 0
            x = 110 + i * 590
            put_in(cv, self.cube(i, q), x, 215, self.c(0, 0.3 + 0.3 * i))
            lab, F, eff = self.labels[i]
            put_in(cv, lab, x + 20, 735, self.c(0, 0.3 + 0.3 * i))
            put_in(cv, F, x + 20, 795, self.c(i + 1))
            put_in(cv, eff, x + 20, 870, self.c(i + 1, 2.0))
        put_in(cv, self.note, 110, 960, self.c(0, 2.0))


class RunVasp(Scene):
    chapter = "Step 3 · Run VASP"
    title = "Step 3 · Run VASP in every folder"
    sentences = [
        ("Copy INCAR, KPOINTS and POTCAR into every strain folder and run VASP.", None, 1.0),
        ("The right settings depend on your material, so do your own convergence tests.", None, 0.3),
        ("Converge the cutoff energy and the k-point mesh until the elastic constants, not just the total energy, stop changing.", None, 0.3),
        ("Use the same settings for every strain, and if the atoms relax inside the strained cell, converge the forces as well.", None, 1.0),
    ]

    def setup(self, D):
        lines = D["incar"].split("\n")
        self.incar_lines = lines
        self.kp = D["kpoints"]
        hl1 = tuple(i for i, line in enumerate(lines) if line.startswith(("ENCUT", "PREC", "EDIFF")))
        hl2 = tuple(i for i, line in enumerate(lines) if line.startswith("ISMEAR"))
        self.incar = [code_panel("INCAR", D["incar"], 820, 26, h=420),
                      code_panel("INCAR", D["incar"], 820, 26, hl1, YELLOW, h=420),
                      code_panel("INCAR", D["incar"], 820, 26, hl2, GREEN, h=420)]
        self.kpts = [code_panel("KPOINTS", self.kp, 820, 26, h=330),
                     code_panel("KPOINTS", self.kp, 820, 26, (3,), YELLOW, h=330)]
        T = self.term = Term(990, 230, 820, 330, size=22)
        t = T.cmd("for d in */strain_*; do", self.cues[0] + 0.3)
        t = T.cmd("  cp INCAR KPOINTS POTCAR $d", t + 0.2, "> ")
        t = T.cmd("  (cd $d && mpirun -np 4 vasp_std)", t + 0.2, "> ")
        T.cmd("done", t + 0.2, "> ")
        card = panel_img(820, 300, fill=(24, 36, 40), border=GREEN)
        d = ImageDraw.Draw(card)
        d.text((30, 24), "Do your own convergence tests", font=sf(30, "Semibold"), fill=GREEN)
        for i, s in enumerate(("cutoff energy (ENCUT)",
                               "k-point mesh and smearing",
                               "converge C, not only E",
                               "same settings for every strain")):
            d.text((34, 84 + i * 52), "•  " + s, font=sf(28), fill=FG)
        self.card = card

    def draw(self, cv, D):
        state = 1 if 0 <= self.c(2) and self.c(3) < 0 else 0
        op, dy = appear(self.t - 0.1)
        put(cv, self.incar[state], 110, 230, op, dy)
        put(cv, self.kpts[1 if 0 <= self.c(2) and self.c(3) < 0 else 0], 110, 680, op, dy)
        self.term.draw(cv, self.t, self.c(0))
        put_in(cv, self.card, 990, 620, self.c(1, 0.2))


class Analyse(Scene):
    chapter = "Step 4 · Analyse"
    title = "Step 4 · Fit the energies"
    sentences = [
        ("When all runs have finished, run elastic-vasp post cubic.", None, 1.6),
        ("It reads the energy from every OUTCAR, checks that each run completed, and fits a parabola to each deformation.", None, 0.6),
        ("Always look at these curves before trusting the numbers: smooth, with the minimum at zero strain.", None, 0.3),
        ("The volume-changing deformation is slightly asymmetric. That is anharmonicity, and the quadratic fit averages it out.", None, 1.0),
    ]

    def setup(self, D):
        T = self.term = Term(110, 215, 1700, 290, size=24)
        t = T.cmd("elastic-vasp post cubic --output C.txt", self.cues[0] + 0.4)
        s = D["std"]
        T.out(["second order coefficients (a in E/V0 = a*delta^2, GPa) are as follow:",
               dict(s="C11_C12_I/=  %.4f" % s["C11_C12_I/"], c=CYAN),
               dict(s="C11_C12_II/= %.4f" % s["C11_C12_II/"], c=ORANGE),
               dict(s="C44/=        %.4f" % s["C44/"], c=PURPLE)], t + 0.6, 0.15)
        self.cache = {}
        self.D = D

    def plot(self, ns, fr, hl):
        key = (ns, tuple(round(f * 12) for f in fr), hl)
        if key not in self.cache:
            fig = plt.figure(figsize=(17, 5.0), dpi=100)
            s = self.D["std"]
            meta = [("C11_C12_I", "a = C11 − C12 = %.2f GPa" % s["C11_C12_I/"], CYAN),
                    ("C11_C12_II", "a = 3/2 C11 + 3 C12 = %.2f GPa" % s["C11_C12_II/"], ORANGE),
                    ("C44", "a = 3/2 C44 = %.2f GPa" % s["C44/"], PURPLE)]
            for i, (name, lab, col) in enumerate(meta):
                ax = fig.add_axes([0.05 + i * 0.33, 0.15, 0.27, 0.68])
                d, e = self.D["curves"][name]
                strain_plot(ax, d, e, ns[i], key[1][i] / 12, col, name)
                ax.text(0.5, 0.9, lab, transform=ax.transAxes, color=hexc(col), fontsize=15, ha="center")
                if hl == i:
                    for sp in ax.spines.values():
                        sp.set_color(hexc(YELLOW))
                        sp.set_linewidth(2.5)
            self.cache[key] = fig_to_sprite(fig)
        return self.cache[key]

    def draw(self, cv, D):
        self.term.draw(cv, self.t, self.t - 0.1)
        ns, fr = [], []
        for i in range(3):
            s0 = self.c(1, 0.3 + 1.7 * i)
            ns.append(int(np.clip(s0 / 0.13 + 1, 0, 7)) if s0 >= 0 else 0)
            f0 = s0 - 1.0
            fr.append(ease(f0 / 0.7) if f0 > 0 else 0.0)
        hl = 1 if self.c(3) >= 0 else None
        put_in(cv, self.plot(tuple(ns), tuple(fr), hl), 60, 530, self.c(1))


class Results(Scene):
    chapter = "Step 5 · Results"
    title = "Step 5 · The elastic constants"
    sentences = [
        ("The three fitted curvatures give three linear equations in C11, C12 and C44.", None, 0.3),
        ("Solving them gives C11 = 106.6, C12 = 63.4, and C44 = 28.9 GPa.",
         "Solving them gives C one one, 106.6, C one two, 63.4, and C four four, 28.9 gigapascals.", 0.4),
        ("All eigenvalues are positive, so the crystal is mechanically stable.", None, 0.3),
        ("The bulk modulus from these constants, 77.8 GPa, is close to the 77.0 GPa from the Birch–Murnaghan fit.", None, 0.3),
        ("They need not agree exactly: the elastic constants are the curvature at one volume and ignore the pressure dependence of the bulk modulus, which the equation of state includes.", None, 1.2),
    ]

    def setup(self, D):
        s = D["std"]
        rows = [("C11_C12_I", r"$C_{11} - C_{12}$", s["C11_C12_I/"], CYAN),
                ("C11_C12_II", r"$\frac{3}{2}C_{11} + 3C_{12}$", s["C11_C12_II/"], ORANGE),
                ("C44", r"$\frac{3}{2}C_{44}$", s["C44/"], PURPLE)]
        tbl = panel_img(820, 330)
        d = ImageDraw.Draw(tbl)
        for j, h in enumerate(("deformation", "a  =", "fitted a (GPa)")):
            d.text((30 + (0, 290, 590)[j], 22), h, font=sf(24, "Semibold"), fill=MUTED)
        d.line((20, 68, 800, 68), fill=BORDER, width=2)
        for i, (name, eq, val, col) in enumerate(rows):
            y = 92 + i * 76
            d.text((30, y + 6), name, font=mono(26, True), fill=col)
            m = math_sprite(eq, 40)
            tbl.alpha_composite(m, (290, y + 22 - m.height // 2 + 6))
            d.text((590, y + 4), "%.2f" % val, font=sf(32, "Medium"), fill=FG)
        self.tbl = tbl
        c11, c12, c44 = s["C11"], s["C12"], s["C44"]
        f = lambda v: "%.1f" % v  # noqa: E731
        M = [[f(c11), f(c12), f(c12), "0", "0", "0"], [f(c12), f(c11), f(c12), "0", "0", "0"],
             [f(c12), f(c12), f(c11), "0", "0", "0"], ["0", "0", "0", f(c44), "0", "0"],
             ["0", "0", "0", "0", f(c44), "0"], ["0", "0", "0", "0", "0", f(c44)]]
        cols = tuple(tuple(MUTED if v == "0" else (CYAN if i == j and i < 3 else PURPLE if i == j else ORANGE)
                           for j, v in enumerate(r)) for i, r in enumerate(M))
        self.mat = matrix_sprite(tuple(tuple(r) for r in M), 30, cols, colw=112)
        self.C = math_sprite(r"$C\;=$", 52)
        self.gpa = text_sprite("GPa", 30, MUTED)
        self.stable = chip("✓  mechanically stable — all eigenvalues positive", 30, GREEN, (20, 44, 36), GREEN, "Semibold")
        K = (c11 + 2 * c12) / 3
        self.check = rich_sprite([("K = (C11 + 2 C12) / 3 = ", FG, "Medium"), ("%.1f GPa" % K, YELLOW, "Bold"),
                                  ("      Birch–Murnaghan: ", FG2, "Regular"), ("%.1f GPa" % D["eos"]["B"], YELLOW, "Bold")], 34)
        self.note = text_sprite("Not identical by construction: the equation of state includes the pressure dependence of B, "
                                "the elastic constants do not.", 25, FG2)
        self.hill = text_sprite("Hill averages:  K = %.2f   G = %.2f   E = %.2f GPa   ν = %.2f" % (
            s["K_VRH"], s["G_VRH"], s["E_VRH"], s["nu_VRH"]), 26, MUTED, is_mono=True)

    def draw(self, cv, D):
        put_in(cv, self.tbl, 110, 240, self.c(0))
        put_in(cv, self.C, 975, 360, self.c(1))
        put_in(cv, self.mat, 1095, 250, self.c(1))
        put_in(cv, self.gpa, 1740, 560, self.c(1, 0.5))
        put_in(cv, self.stable, 110, 640, self.c(2, 0.2))
        put_in(cv, self.check, 110, 760, self.c(3, 0.2))
        put_in(cv, self.note, 110, 822, self.c(4, 0.2))
        put_in(cv, self.hill, 110, 885, self.c(3, 2.5))


class Auto(Scene):
    chapter = "Any orientation"
    title = "Any orientation: the automatic mode"
    sentences = [
        ("What if your cell is not in the standard orientation?", None, 0.4),
        ("We rotated the same cell by arbitrary angles and ran the automatic mode: auto pre, then auto post.", None, 0.6),
        ("It found the cubic symmetry, chose three well-conditioned deformations, and rotated the result back to the standard setting.", None, 0.3),
        ("Both routes agree within the numerical noise of about one to one and a half GPa.", None, 1.2),
    ]

    def setup(self, D):
        self.cells = [self.cell(D["lat_std"], CYAN, "standard cell"), self.cell(D["lat_rot"], ORANGE, "rotated 25°, 40°, 10°")]
        T = self.term = Term(960, 230, 850, 400, size=21)
        t = T.cmd("elastic-vasp pre auto --max-strain 0.03 --points 7", self.cues[1] + 0.5)
        t = T.out([dict(s="Crystal symmetry (spglib, symprec=0.001 A): Fm-3m (225)", hl=(self.cues[2] + 0.2, CYAN)),
                   "Deformations (3):",
                   dict(s="    e1+e3    (1, 0, 1, 0, 0, 0) * delta", hl=(self.cues[2] + 1.6, ORANGE)),
                   dict(s="    e2-e3    (0, 1, -1, 0, 0, 0) * delta", hl=(self.cues[2] + 1.8, ORANGE)),
                   dict(s="    e1+e4    (1, 0, 0, 1, 0, 0) * delta", hl=(self.cues[2] + 2.0, ORANGE))], t + 0.3, 0.1)
        T.cmd("elastic-vasp post auto --frame standard", t + 0.6)
        r, s = D["rot"], D["std"]
        tbl = panel_img(1700, 250)
        d = ImageDraw.Draw(tbl)
        for j, h in enumerate(("C11", "C12", "C44")):
            d.text((1000 + j * 230, 22), h, font=sf(30, "Semibold"), fill=MUTED)
        d.line((20, 74, 1680, 74), fill=BORDER, width=2)
        for i, (lab, col, vals) in enumerate((("cubic_pre / cubic_post — standard cell", CYAN, (s["C11"], s["C12"], s["C44"])),
                                              ("auto_pre / auto_post — rotated cell", ORANGE, (r["C11"], r["C12"], r["C44"])))):
            y = 98 + i * 72
            d.text((30, y), lab, font=sf(32, "Medium"), fill=col)
            for j, v in enumerate(vals):
                d.text((1000 + j * 230, y), "%.1f" % v, font=sf(34, "Semibold"), fill=FG)
        d.text((1690 - sf(24).getlength("GPa") - 10, 22), "GPa", font=sf(24), fill=MUTED)
        self.tbl = tbl

    def cell(self, lat, col, label):
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
        fig = plt.figure(figsize=(4.0, 4.3), dpi=100)
        ax = fig.add_axes([0, 0.08, 1, 0.92], projection="3d")
        ax.set_axis_off()
        ax.patch.set_alpha(0)
        a1, a2, a3 = lat
        corners = [i * a1 + j * a2 + k * a3 for i in (0, 1) for j in (0, 1) for k in (0, 1)]
        c = np.array(corners)
        for i in range(8):
            for j in range(i + 1, 8):
                if bin(i ^ j).count("1") == 1:
                    ax.plot(*zip(c[i], c[j]), color=hexc(col), lw=2.4)
        ax.scatter(c[:, 0], c[:, 1], c[:, 2], s=70, color=hexc(FG2), depthshade=False)
        for v, name in ((np.eye(3)[0], "x"), (np.eye(3)[1], "y"), (np.eye(3)[2], "z")):
            ax.quiver(-1.5, -1.5, -1.5, *(v * 1.6), color=hexc(MUTED), lw=1.5, arrow_length_ratio=0.25)
            ax.text(*(np.array([-1.5] * 3) + v * 2.0), name, color=hexc(MUTED), fontsize=13)
        ctr = c.mean(0)
        ax.set_xlim(ctr[0] - 3, ctr[0] + 3)
        ax.set_ylim(ctr[1] - 3, ctr[1] + 3)
        ax.set_zlim(ctr[2] - 3, ctr[2] + 3)
        ax.view_init(elev=18, azim=-60)
        fig.text(0.5, 0.02, label, ha="center", color=hexc(col), fontsize=16)
        return fig_to_sprite(fig)

    def draw(self, cv, D):
        put_in(cv, self.cells[0], 90, 210, self.c(0))
        put_in(cv, self.cells[1], 500, 210, self.c(1))
        self.term.draw(cv, self.t, self.c(1))
        put_in(cv, self.tbl, 110, 700, self.c(3))


class Visualise(Scene):
    chapter = "Visualise"
    title = "Visualise the anisotropy"
    sentences = [
        ("To see how the stiffness depends on direction, install the companion package elastic-vis from GitHub.", None, 0.6),
        ("Save the matrix with --output, and pass it to elastic-vis.", None, 0.4),
        ("For aluminium, Young's modulus ranges from 59.4 GPa along ⟨100⟩ to 77.2 GPa along ⟨111⟩.",
         "For aluminium, Young's modulus ranges from 59.4 gigapascals along the cube axes, to 77.2 along the body diagonal.", 1.4),
    ]

    def setup(self, D):
        T = self.term = Term(110, 215, 1700, 190, size=24)
        T.cmd('pip install "git+https://github.com/prabhat7iitm/elastic-vis.git"', self.cues[0] + 0.5)
        T.cmd("elastic-vis C.txt --plot --Young --Shear --2D --3D --format png --output-prefix Al", self.cues[1] + 0.3)
        base = os.path.join(REPO, "docs/_static/elastic_vis")
        self.i3d = image_card(os.path.join(base, "Al_Youngs_Modulus_3D.png"), 560, 580)
        self.i2d = image_card(os.path.join(base, "Al_Youngs_Modulus_Polar_2D.png"), 1100, 430)
        self.e1 = chip("E<100> = 59.4 GPa", 32, CYAN, is_mono=True)
        self.e2 = chip("E<111> = 77.2 GPa", 32, ORANGE, is_mono=True)

    def draw(self, cv, D):
        self.term.draw(cv, self.t, self.t - 0.1)
        put_in(cv, self.i3d, 110, 440, self.c(1, 3.0))
        x = 110 + self.i3d.width + 40
        put_in(cv, self.i2d, x, 440, self.c(1, 3.7))
        y = 440 + self.i2d.height + 45
        put_in(cv, self.e1, x, y, self.c(2, 0.8))
        put_in(cv, self.e2, x + self.e1.width + 40, y, self.c(2, 2.6))


class Lessons(Scene):
    chapter = "Lessons"
    title = "Lessons from real runs"
    sentences = [
        ("A few lessons from getting here.", None, 0.2),
        ("A Γ-point-only calculation of a metal gives nonsense: here the energy even went down under strain.", None, 0.3),
        ("An unconverged k-point mesh gave noise as large as the signal: that is why the convergence tests matter.", None, 0.3),
        ("So do your own convergence tests, and always look at the energy–strain curves.", None, 1.2),
    ]

    def setup(self, D):
        def card(title, mark, col, lines, bg):
            img = panel_img(540, 470, fill=bg, border=col)
            d = ImageDraw.Draw(img)
            d.text((30, 26), mark, font=sf(54, "Bold"), fill=col)
            d.text((100, 36), title, font=sf(32, "Semibold"), fill=FG)
            y = 130
            for line in lines:
                for k, w in enumerate(wrap(line, sf(28), 470)):
                    d.text((40, y), ("•  " if k == 0 else "    ") + w, font=sf(28), fill=FG2)
                    y += 40
                y += 16
            return img
        self.cards = [
            card("Γ point only", "✗", RED, ["Metal: degenerate states at the Fermi level", "Energy decreased under strain",
                                           "Constants meaningless"], (40, 22, 32)),
            card("Unconverged k-mesh", "✗", RED, ["Mesh too coarse for a metal", "~1 meV k-point noise",
                                                 "As large as the signal at ±2 %"], (40, 22, 32)),
            card("Do your own tests", "✓", GREEN, ["Run your own convergence tests", "Look at the curves!"], (20, 40, 34)),
        ]

    def draw(self, cv, D):
        for i, cd in enumerate(self.cards):
            put_in(cv, cd, 110 + i * 590, 300, self.c(i + 1))


class Outro(Scene):
    chapter = "Learn more"
    title = ""
    sentences = [
        ("The full documentation, with this worked example, is on the project website.", None, 0.3),
        ("If elastic_vasp helps your work, please cite the papers listed there. Thanks for watching!", None, 2.5),
    ]

    def setup(self, D):
        self.t1 = text_sprite("Learn more", 76, FG, "Bold")
        self.u1 = chip("prnvrvs.github.io/elastic_vasp", 40, CYAN, PANEL, CYAN, "Semibold", is_mono=True, pad=(30, 16))
        self.u2 = chip("github.com/prnvrvs/elastic_vasp", 34, FG2, PANEL, BORDER, "Regular", is_mono=True, pad=(26, 12))
        self.cite = text_sprite(
            "Please cite:\n"
            "P. Kumar and I. Adlakha, J. Eng. Mater. Technol. 145 (1), 2022 — doi:10.1115/1.4055097\n"
            "P. Mishra, P. Kumar, L. Neelakantan and I. Adlakha, Comput. Mater. Sci. 214, 111667 (2022)",
            27, FG2, "Regular", spacing=1.5)
        self.thanks = text_sprite("Thanks for watching!", 46, YELLOW, "Semibold")
        self.img = image_card(os.path.join(REPO, "docs/_static/elastic_vis/Al_Youngs_Modulus_3D.png"), 560, 560)

    def draw(self, cv, D):
        put_in(cv, self.t1, 130, 200, self.t - 0.1)
        put_in(cv, self.u1, 130, 340, self.c(0, 0.2))
        put_in(cv, self.u2, 130, 450, self.c(0, 1.0))
        put_in(cv, self.cite, 130, 600, self.c(1, 0.2))
        put_in(cv, self.thanks, 130, 850, self.c(1, 3.5))
        put_in(cv, self.img, 1250, 180, self.t - 0.4, 0.8)


SCENES = [Intro, Idea, Install, Structure, Check, Generate, Deform, RunVasp, Analyse, Results, Auto, Visualise, Lessons, Outro]


# --------------------------------------------------------------------------- frame furniture

def background():
    y = np.linspace(0, 1, H)[:, None]
    x = np.linspace(0, 1, W)[None, :]
    base = np.array(BG1)[None, None, :] * (1 - y[..., None]) + np.array(BG2)[None, None, :] * y[..., None]
    glow = np.exp(-(((x - 0.85) / 0.35) ** 2 + ((y - 0.1) / 0.35) ** 2))[..., None] * np.array([20, 40, 70])
    arr = np.clip(base + glow, 0, 255).astype(np.uint8)
    return Image.fromarray(arr, "RGB").convert("RGBA")


def scene_base(bg, scene, idx, n):
    img = bg.copy()
    d = ImageDraw.Draw(img)
    d.text((110, 52), "elastic_vasp", font=mono(24, True), fill=CYAN)
    d.text((110 + mono(24, True).getlength("elastic_vasp") + 14, 54), "·  tutorial", font=sf(23), fill=MUTED)
    lab = "%02d / %02d   %s" % (idx + 1, n, scene.chapter)
    d.text((W - 110 - sf(23).getlength(lab), 54), lab, font=sf(23), fill=MUTED)
    if scene.title:
        d.text((108, 112), scene.title, font=sf(54, "Bold"), fill=FG)
    return img


def progress(cv, t, total, starts):
    d = ImageDraw.Draw(cv)
    y = H - 34
    d.rectangle((110, y, W - 110, y + 4), fill=(40, 50, 78))
    d.rectangle((110, y, 110 + (W - 220) * t / total, y + 4), fill=CYAN)
    for s in starts[1:]:
        x = 110 + (W - 220) * s / total
        d.rectangle((x - 1, y - 3, x + 1, y + 7), fill=(90, 104, 140))


# --------------------------------------------------------------------------- subtitles and description

def fmt_srt(t):
    ms = int(round(t * 1000))
    return "%02d:%02d:%02d,%03d" % (ms // 3600000, ms // 60000 % 60, ms // 1000 % 60, ms % 1000)


def srt_chunks(text, start, dur, maxc=84):
    words, chunks, cur = text.split(), [], ""
    maxc = len(text) / -(-len(text) // maxc) + 6
    for w in words:
        if len(cur) + len(w) + 1 > maxc and cur:
            chunks.append(cur)
            cur = w
        else:
            cur = (cur + " " + w).strip()
    chunks.append(cur)
    total = sum(len(c) for c in chunks)
    out, t = [], start
    for c in chunks:
        d = dur * len(c) / total
        lines = wrap(c, sf(30), 30 * 0.5 * 44)
        out.append((t, t + d, "\n".join(lines)))
        t += d
    return out


def mmss(t):
    return "%d:%02d" % (int(t) // 60, int(t) % 60)


# --------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "build"))
    ap.add_argument("--engine", choices=("kokoro", "chatterbox", "say"), default="kokoro")
    ap.add_argument("--voice-sample", help="chatterbox: 10-20 s recording of the voice to clone")
    ap.add_argument("--exaggeration", type=float, default=0.5, help="chatterbox: expressiveness")
    ap.add_argument("--cfg", type=float, default=0.5, help="chatterbox: pacing / guidance weight")
    ap.add_argument("--voice", help="Kokoro voice (default am_michael) or macOS voice (default Daniel)")
    ap.add_argument("--speed", type=float, default=1.0, help="Kokoro speaking speed")
    ap.add_argument("--rate", type=int, default=172, help="macOS say words per minute")
    ap.add_argument("--tts-python", default=sys.executable, help="Python with kokoro-onnx installed")
    ap.add_argument("--tts-models", help="folder with kokoro-v1.0.onnx and voices-v1.0.bin")
    ap.add_argument("--recorded", help="folder with your own recordings s01_1.m4a ... (ids in narration.txt)")
    ap.add_argument("--ffmpeg", default="ffmpeg", help="ffmpeg executable (needs libx264)")
    ap.add_argument("--only", type=int, nargs="*", help="render only these scene numbers (1-based) as stills")
    args = ap.parse_args()
    args.voice = args.voice or ("am_michael" if args.engine == "kokoro" else "Daniel")
    if args.engine == "kokoro" and not args.tts_models:
        ap.error("--tts-models is required with the kokoro engine")
    if args.engine == "chatterbox":
        if not args.voice_sample:
            ap.error("--voice-sample is required with the chatterbox engine")
        h = zlib.crc32(open(args.voice_sample, "rb").read())
        args.voice = "clone%08x" % h
    os.makedirs(os.path.join(args.out, "voice"), exist_ok=True)
    D = read_data(args.data)

    # narration: one clip per sentence, exact cue times
    scenes = [cls() for cls in SCENES]
    apply_script(scenes)
    pcm = bytearray()
    rate = 24000
    silence = lambda s: b"\x00\x00" * int(round(s * rate))  # noqa: E731
    tag = "%s_%s" % (args.voice.split()[0], args.speed if args.engine == "kokoro" else args.rate
                     if args.engine == "say" else "%g_%g" % (args.exaggeration, args.cfg))
    clip = lambda k, j: os.path.join(args.out, "voice", "s%02d_%d_%s_%08x.wav" % (  # noqa: E731
        k, j, tag, zlib.crc32(spoken(scenes[k].sentences[j][1] or scenes[k].sentences[j][0]).encode())))
    synth_all([(spoken(say_text or sub), clip(k, j)) for k, sc in enumerate(scenes)
               for j, (sub, say_text, pause) in enumerate(sc.sentences)], args)
    if args.recorded:
        tts_clip, clip = clip, lambda k, j: recorded_clip(args, k, j) or tts_clip(k, j)  # noqa: E731
    subs, t_global, starts = [], 0.0, []
    for k, sc in enumerate(scenes):
        starts.append(t_global)
        sc.cues = []
        t = LEAD
        pcm += silence(LEAD)
        for j, (sub, say_text, pause) in enumerate(sc.sentences):
            frames, dur = read_clip(clip(k, j))
            sc.cues.append(t)
            subs += srt_chunks(sub, t_global + t, dur)
            pcm += frames
            t += dur
            gap = GAP + pause
            pcm += silence(gap)
            t += gap
        sc.dur = t + TAIL - GAP
        pcm += silence(TAIL - GAP)
        n_samples = int(round((t_global + sc.dur) * rate))
        pcm = pcm[:2 * n_samples] + silence(max(0, n_samples - len(pcm) // 2) / rate)
        t_global += sc.dur
        sc.setup(D)
    total = t_global

    with wave.open(os.path.join(args.out, "narration.wav"), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(rate)
        w.writeframes(bytes(pcm))
    with open(os.path.join(args.out, "tutorial.srt"), "w") as f:
        for i, (a, b, text) in enumerate(subs, 1):
            f.write("%d\n%s --> %s\n%s\n\n" % (i, fmt_srt(a), fmt_srt(b - 0.05), text))
    write_description(args.out, scenes, starts, total, D)

    bg = background()
    bases = [scene_base(bg, sc, i, len(scenes)) for i, sc in enumerate(scenes)]

    def frame(k, t):
        sc = scenes[k]
        sc.t = t
        cv = bases[k].copy()
        sc.draw(cv, D)
        progress(cv, starts[k] + t, total, starts)
        if t < 0.4:
            cv = Image.blend(bg, cv, ease(t / 0.4))
        return cv

    if args.only:
        for k in args.only:
            sc = scenes[k - 1]
            for frac in (0.35, 0.999):
                frame(k - 1, sc.dur * frac).convert("RGB").save(os.path.join(args.out, "still_%02d_%d.png" % (k, int(frac * 100))))
        return

    thumbnail(args.out, D)
    cmd = [args.ffmpeg, "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", "%dx%d" % (W, H),
           "-r", str(FPS), "-i", "-", "-i", os.path.join(args.out, "narration.wav"),
           "-c:v", "libx264", "-preset", "medium", "-crf", "18", "-pix_fmt", "yuv420p", "-tune", "stillimage",
           "-af", "loudnorm=I=-16:TP=-1.5:LRA=11", "-c:a", "aac", "-b:a", "160k", "-ar", "48000", "-shortest", "-movflags", "+faststart",
           os.path.join(args.out, "tutorial.mp4")]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    for k, sc in enumerate(scenes):
        nf = int(round(sc.dur * FPS))
        for i in range(nf):
            proc.stdin.write(frame(k, i / FPS).convert("RGB").tobytes())
        print("scene %2d/%d done (%.1f s)" % (k + 1, len(scenes), sc.dur), flush=True)
    proc.stdin.close()
    if proc.wait() != 0:
        raise SystemExit("ffmpeg failed (exit code %d)" % proc.returncode)
    print("total %.1f s -> %s" % (total, os.path.join(args.out, "tutorial.mp4")))


def write_description(out, scenes, starts, total, D):
    s = D["std"]
    chapters = "\n".join("%s %s" % (mmss(t), sc.chapter) for t, sc in zip(starts, scenes))
    text = """Elastic constants with VASP: a step-by-step tutorial for elastic_vasp (fcc Al, real data)

How to compute the second-order elastic constants of a crystal from VASP energy–strain calculations with the open-source Python package elastic_vasp. Every number in the video comes from real VASP (PBE) calculations on fcc aluminium.

Result: C11 = {c11:.1f}, C12 = {c12:.1f}, C44 = {c44:.1f} GPa; bulk modulus {k:.1f} GPa (Birch–Murnaghan equation of state: {b:.1f} GPa).

Chapters
{chapters}

Documentation: https://prnvrvs.github.io/elastic_vasp/
Worked example: https://prnvrvs.github.io/elastic_vasp/worked_example.html
Source code: https://github.com/prnvrvs/elastic_vasp
Install: pip install "git+https://github.com/prnvrvs/elastic_vasp.git"
Visualisation (elastic-vis): pip install "git+https://github.com/prabhat7iitm/elastic-vis.git"

Please cite
1. P. Kumar and I. Adlakha, "Effect of Interstitial Hydrogen on Elastic Behavior of Metals: An Ab-Initio Study", J. Eng. Mater. Technol. 145 (1), 2022. https://doi.org/10.1115/1.4055097
2. P. Mishra, P. Kumar, L. Neelakantan and I. Adlakha, "First-principles prediction of electrochemical polarization and mechanical behavior in Mg based intermetallics", Comput. Mater. Sci. 214, 111667 (2022). https://doi.org/10.1016/j.commatsci.2022.111667

#VASP #DFT #ElasticConstants #ComputationalMaterialsScience #Python
""".format(c11=s["C11"], c12=s["C12"], c44=s["C44"], k=(s["C11"] + 2 * s["C12"]) / 3, b=D["eos"]["B"], chapters=chapters)
    with open(os.path.join(out, "youtube_description.txt"), "w") as f:
        f.write(text)


def thumbnail(out, D):
    bg = background()
    img = bg.copy()
    d = ImageDraw.Draw(img)
    d.text((110, 150), "ELASTIC", font=sf(170, "Heavy"), fill=FG)
    d.text((110, 330), "CONSTANTS", font=sf(170, "Heavy"), fill=FG)
    d.text((116, 530), "with VASP", font=sf(110, "Bold"), fill=CYAN)
    d.text((118, 700), "elastic_vasp  ·  step by step  ·  real data", font=sf(48, "Medium"), fill=FG2)
    s = D["std"]
    fig = plt.figure(figsize=(7.6, 6.2), dpi=100)
    ax = fig.add_axes([0.12, 0.12, 0.84, 0.84])
    for name, col in (("C11_C12_II", ORANGE), ("C11_C12_I", CYAN), ("C44", PURPLE)):
        dd, e = D["curves"][name]
        style_ax(ax)
        p = np.polyfit(dd, e, 2)
        x = np.linspace(-3.3, 3.3, 200)
        ax.plot(x, np.polyval(p, x), color=hexc(col), lw=4)
        ax.plot(dd, e, "o", ms=13, color=hexc(FG), mec=hexc(col), mew=3)
    ax.set_xlabel("strain δ (%)", fontsize=20)
    ax.set_ylabel("E − E₀ (meV)", fontsize=20)
    ax.tick_params(labelsize=16)
    img.alpha_composite(fig_to_sprite(fig), (1140, 120))
    badge = chip("C11 = %.1f   C12 = %.1f   C44 = %.1f GPa" % (s["C11"], s["C12"], s["C44"]), 44, YELLOW,
                 (45, 40, 20), YELLOW, "Bold", pad=(34, 18))
    img.alpha_composite(badge, (110, 860))
    img.convert("RGB").resize((1280, 720), Image.LANCZOS).save(os.path.join(out, "thumbnail.png"))


if __name__ == "__main__":
    main()
