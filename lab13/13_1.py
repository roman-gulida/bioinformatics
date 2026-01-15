import numpy as np


def predict_n_states(matrix, initial_vector, steps=5):
    predictions = [initial_vector]
    current_state = initial_vector

    for _ in range(steps):
        current_state = matrix @ current_state
        predictions.append(current_state)

    return predictions


def main():
    matrix = np.array([[0.7, 0.2, 0.1], [0.3, 0.5, 0.2], [0.1, 0.3, 0.6]])

    initial_vector = np.array([1.0, 0.0, 0.0])

    print("Transition Matrix:")
    print(matrix)
    print("\nInitial Vector:")
    print(initial_vector)

    predictions = predict_n_states(matrix, initial_vector, steps=5)

    for i, state in enumerate(predictions):
        print(f"\nStep {i}:")
        print(state)
        print(f"Sum: {np.sum(state):.6f}")


if __name__ == "__main__":
    main()
