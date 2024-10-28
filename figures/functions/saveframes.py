from matplotlib.animation import FileMovieWriter
class BunchOFiles(FileMovieWriter):
    supported_formats = ['png', 'jpeg', 'bmp', 'svg', 'pdf']

    def __init__(self, *args, extra_args=None, **kwargs):
        # extra_args aren't used but we need to stop None from being passed
        super().__init__(*args, extra_args=(), **kwargs)

    def setup(self, fig, dpi, frame_prefix):
        super().setup(fig, dpi, frame_prefix)
        self.fname_format_str = '%s%%d.%s'
        self.temp_prefix, self.frame_format = self.outfile.split('.')

    def grab_frame(self, **savefig_kwargs):
        '''
        Grab the image information from the figure and save as a movie frame.
        All keyword arguments in savefig_kwargs are passed on to the 'savefig'
        command that saves the figure.
        '''

        # Tell the figure to save its data to the sink, using the
        # frame format and dpi.
        with self._frame_sink() as myframesink:
            self.fig.savefig(myframesink, format=self.frame_format,
                              dpi=self.dpi, **savefig_kwargs)

    def finish(self):
        self._frame_sink().close()

