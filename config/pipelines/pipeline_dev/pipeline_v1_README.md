# Pipeline v1: Data Processing and Model Training Pipeline

This pipeline is designed to process data, prepare it for machine learning (ML), and run model training and testing. The pipeline includes several key stages, from processing settings to running and saving a trained machine learning model.

## Steps in the Pipeline

1. **Process Settings**
   - The pipeline first converts the YAML settings file (`settings.yaml`) into a Python format (`processed_settings.py`), which can be easily loaded by the pipeline.

2. **Load Settings**
   - The processed settings are then loaded into Python variables from the `processed_settings.py` file.

3. **Fetch and Preprocess Data**
   - The pipeline fetches data based on the loaded settings and performs basic preprocessing, such as:
     - Preprocessing protein sequences.
     - Preprocessing activity values.
   - The preprocessed data is saved into a CSV file (`processed_data.csv`).

4. **Prepare Data for ML**
   - The fetched and preprocessed data is then prepared for machine learning, transforming it into the right format for input into the ML model. This step generates the `ml_input.csv` file.

5. **Train and Test Model**
   - The ML model is trained and tested using the prepared data (`ml_input.csv`).
   - A model accuracy threshold is defined (e.g., 85%), and the model is saved if the accuracy is achieved.
   - The trained model is stored as `trained_model.pkl`.

6. **Save Model and Results**
   - Once the model has been trained successfully, the pipeline saves the model, accuracy results, and plots (e.g., accuracy plot, model summary).

## How to Use

1. **Configuration Files**:
   - The pipeline requires several configuration files. Ensure that the `settings.yaml` file exists and contains the necessary settings for the pipeline.
   - The settings will be processed and loaded into `processed_settings.py` before running the pipeline.

2. **Running the Pipeline**:
   - You can run the pipeline manually. For example, if you have a script that orchestrates the pipeline, you can call each step sequentially using the settings and data files.

3. **Manual Trigger**:
   - The pipeline is triggered manually using the command line. After setting up the environment and necessary configuration files, execute the pipeline script to run all the steps.

4. **Model Training**:
   - In the "Train and Test Model" step, the model will be trained and tested. If the accuracy meets the threshold, the model will be saved to `models/trained_model.pkl`.
   - The results (accuracy, plots, model summary) will be stored in the `models/` directory.

## Files Created

- `data/processed_data.csv`: Preprocessed data ready for machine learning.
- `data/ml_input.csv`: Data formatted for ML model input.
- `models/trained_model.pkl`: The trained model.
- `models/model_results.json`: Results from the model training (e.g., accuracy, loss).
- `models/accuracy_plot.png`: Accuracy plot.
- `models/model_summary.txt`: Text file with a summary of the trained model.

## Requirements

- Python 3.x
- Required Python packages (see `requirements.txt`)

## Execution Example

To run the pipeline manually, you might run something like:

```bash
python src/pipeline/run_pipeline.py --config=config/settings.yaml
