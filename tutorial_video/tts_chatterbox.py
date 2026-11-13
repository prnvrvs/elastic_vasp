"""Batch voice-cloned text-to-speech with Chatterbox (run by make_video.py in a Python env
that has chatterbox-tts).

    python tts_chatterbox.py REFERENCE.wav EXAGGERATION CFG JOBS.json

REFERENCE is a 10-20 s recording of the voice to copy. JOBS.json is a list of
[text, wav_path]; each clip is written as 16-bit mono 24 kHz WAVE.
"""
import json
import os
import sys

import torch
import torchaudio
from chatterbox.tts import ChatterboxTTS

reference, exaggeration, cfg, jobs = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
device = "mps" if torch.backends.mps.is_available() else "cpu"
_load = torch.load
torch.load = lambda *a, **k: _load(*a, **{**k, "map_location": torch.device(device)})  # checkpoints saved on CUDA
model = ChatterboxTTS.from_pretrained(device=device)
torch.manual_seed(0)
for text, path in json.load(open(jobs)):
    wav = model.generate(text, audio_prompt_path=reference, exaggeration=exaggeration, cfg_weight=cfg)
    wav = wav.detach().cpu().float()
    if model.sr != 24000:
        wav = torchaudio.functional.resample(wav, model.sr, 24000)
    wav = wav.reshape(1, -1).clamp(-1, 1)
    torchaudio.save(path, wav, 24000, encoding="PCM_S", bits_per_sample=16)
    print(os.path.basename(path), "%.1f s" % (wav.shape[1] / 24000), flush=True)
