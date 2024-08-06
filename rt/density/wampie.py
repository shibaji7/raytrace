class WAMIPE2d(object):
    def __init__(
        self,
        cfg,
        event,
    ):
        self.cfg = cfg
        self.file_name = self.cfg.density_file_location
        self.event = event
        return
