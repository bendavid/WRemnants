import os, logging

class CustomFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""

    green = "\x1b[1;32m"
    grey = "\x1b[38;21m"
    yellow = "\x1b[33;21m"
    red = "\x1b[31;21m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    #format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"
    myformat = "%(levelname)s:%(filename)s: %(message)s"

    FORMATS = {
        logging.DEBUG: green + myformat + reset,
        logging.INFO: grey + myformat + reset,
        logging.WARNING: yellow + myformat + reset,
        logging.ERROR: red + myformat + reset,
        logging.CRITICAL: bold_red + myformat + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

logging_verboseLevel = [logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG]

def set_logging_level(log, verbosity):
    log.setLevel(logging_verboseLevel[max(0, min(4, verbosity))])

def setup_logger(basefile, verbosity=3, no_colors=False):
    setup_func = setup_base_logger if no_colors else setup_color_logger
    return setup_func(os.path.basename(basefile), verbosity)

def setup_color_logger(name, verbosity):
    base_logger = logging.getLogger("wremnants")
    # set console handler
    ch = logging.StreamHandler()
    ch.setFormatter(CustomFormatter())
    base_logger.addHandler(ch)
    set_logging_level(base_logger, verbosity)
    base_logger.propagate = False # to avoid propagating back to root logger, which would print messages twice
    return base_logger.getChild(name)
    
def setup_base_logger(name, verbosity):
    logging.basicConfig(format='%(levelname)s: %(message)s')
    base_logger = logging.getLogger("wremnants")
    set_logging_level(base_logger, verbosity)
    return base_logger.getChild(name)
    
def child_logger(name):
    return logging.getLogger("wremnants").getChild(name)
