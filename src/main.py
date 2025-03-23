from pipeline.fetch_data import fetch_data
from pipeline.process_data import process_data
from pipeline.store_data import store_data
from models.train_model import train_model
from utils.logger import setup_logger

def main():
    # Set up logging

    # Set up the logger (you can adjust the log level or specify a file to save logs)
    logger = setup_logger(name="my_application", log_level=logging.DEBUG, log_file="app.log")

    # Example logging messages
    logger.debug("This is a debug message")
    logger.info("This is an informational message")
    logger.warning("This is a warning message")
    logger.error("This is an error message")
    logger.critical("This is a critical message")


    # Run data pipeline
    data = fetch_data()
    processed_data = process_data(data)
    store_data(processed_data)

    # Train model
    model = train_model(processed_data)

    # Log success
    logger.info("Pipeline and model training completed successfully.")

if __name__ == "__main__":
    main()
