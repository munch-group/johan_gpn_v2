def prediction_step(self, model, inputs, prediction_loss_only, ignore_keys=None):
    # Debugging: Check the inputs
    print("Debug: Checking inputs in prediction_step")
    print(f"Inputs type: {type(inputs)}")
    print(f"Inputs content: {inputs}")

    if inputs is None:
        raise ValueError("The 'inputs' variable is None. Ensure the dataset is properly prepared.")

    # ...existing code...