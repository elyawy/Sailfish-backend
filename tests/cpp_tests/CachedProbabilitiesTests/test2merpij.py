import numpy as np
from itertools import product

AMINO_ACIDS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

data = np.loadtxt('tests/cpp_tests/CachedProbabilitiesTests/somePij.csv', delimiter=',')


mers_2 = list("".join(i) for i in product(AMINO_ACIDS, repeat=2))


probabilities = np.zeros((len(mers_2), len(mers_2)))


sum_probs = 0.0

for first in mers_2:
    for second in mers_2:
        parent_char_0 = AMINO_ACIDS.index(first[0])
        next_char_0 = AMINO_ACIDS.index(second[0])
        parent_char_1 = AMINO_ACIDS.index(first[1])
        next_char_1 = AMINO_ACIDS.index(second[1])
        probabilities[mers_2.index(first), mers_2.index(second)] = data[parent_char_0, next_char_0] * data[parent_char_1, next_char_1]
        cpp_value = data[parent_char_0, next_char_0] * data[parent_char_1, next_char_1]
        # Here you would compute the expected value using your Python model
        # For demonstration, let's assume we have a function `compute_pij`
        # expected_value = compute_pij(first, second, branch_length, rate_category)
        # For now, we'll just print the indices and the C++ value
        sum_probs += cpp_value
        # print(f"Pij({first}, {second}) = {cpp_value}")

print(f"Sum of probabilities: {sum_probs}")

probs_array = np.array(probabilities)

#caluculate the sum of each row for the first 230 columns, and sort by highest sum to lowest
row_sums = probs_array[:, :230].sum(axis=1)
sorted_indices = np.argsort(-row_sums)
sorted_probs = probs_array[sorted_indices, :230]

print("Top 10 rows with highest sums:")
for i in range(230):
    print(f"Row {sorted_indices[i]} sum: {row_sums[sorted_indices[i]]}")


#create heatmap of probabilities with seaborn with log scale
# import seaborn as sns
# import matplotlib.pyplot as plt
# plt.figure(figsize=(12, 10))
# sns.heatmap(np.log(probabilities), xticklabels=mers_2, yticklabels=mers_2, cmap="YlGnBu")
# plt.title("2-mer Transition Probabilities Heatmap")
# plt.xlabel("To 2-mer")
# plt.ylabel("From 2-mer")
# plt.show()
# print(data)
# print(data.shape)