import logging

def setup_logger(name="main_logger", log_level=logging.INFO, log_file=None):
    """
    Set up and return a logger.

    :param name: Name of the logger.
    :param log_level: Log level (default is INFO).
    :param log_file: Log file path (if None, logs will be output to the console).
    :return: Configured logger instance.
    """
    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(log_level)

    # Create formatter to control log message format
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(log_format)

    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # If a log file is specified, create a file handler
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger
