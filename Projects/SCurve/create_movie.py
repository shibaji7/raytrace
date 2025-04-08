import ffmpeg

folder = "/home/shibaji/OneDrive/trace/outputs/April2024_SAMI3_eclipse_hamsci_05MHz_SCurve/2024-04-08/wwv/sami3/w2naf/"
(
    ffmpeg.input(folder + "*.png", pattern_type="glob", framerate=25)
    .output(
        "movie.mp4", crf=20, preset="slower", movflags="faststart", pix_fmt="yuv420p"
    )
    .run()
)
