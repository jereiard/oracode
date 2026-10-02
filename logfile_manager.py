class LogfileManager:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(LogfileManager, cls).__new__(cls)
            cls._logfile = None
        return cls._instance

    @property
    def logfile(self):
        return self._logfile

    @logfile.setter
    def logfile(self, path):
        self._logfile = path

    def close(self):
        # Compatibility shim: reporters open and close file handles per write.
        return None

logfile_manager = LogfileManager()
