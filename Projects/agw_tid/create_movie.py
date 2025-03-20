import ffmpeg

run_name = "May2017_gemini_tid_cosmic2"
out_fname = f"figures/{run_name}.mp4"
folder = f"/home/shibaji/OneDrive/trace/outputs/{run_name}/2017-05-27/fhe/11/gemini/"
(
    ffmpeg
    .input(folder+"*.png", pattern_type="glob", framerate=25)
    .output(out_fname, crf=20, preset="slower", movflags="faststart", pix_fmt="yuv420p")
    .run()
)