"""Batch text-to-speech with Kokoro (run by make_video.py in a Python env that has kokoro-onnx).

    python tts_kokoro.py MODEL_DIR VOICE SPEED JOBS.json

JOBS.json is a list of [text, wav_path]; each clip is written as 16-bit mono 24 kHz WAVE.
MODEL_DIR holds kokoro-v1.0.onnx and voices-v1.0.bin
(https://github.com/thewh1teagle/kokoro-onnx/releases/tag/model-files-v1.0).
"""
import json
import os
import sys

import soundfile as sf
from kokoro_onnx import Kokoro

model_dir, voice, speed, jobs = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
kokoro = Kokoro(os.path.join(model_dir, "kokoro-v1.0.onnx"), os.path.join(model_dir, "voices-v1.0.bin"))
lang = "en-gb" if voice.startswith("b") else "en-us"
for text, path in json.load(open(jobs)):
    audio, sr = kokoro.create(text, voice=voice, speed=speed, lang=lang)
    assert sr == 24000, sr
    sf.write(path, audio, sr, subtype="PCM_16")
    print(os.path.basename(path), "%.1f s" % (len(audio) / sr), flush=True)
