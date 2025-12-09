import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from matplotlib.lines import Line2D

y_axis = ['ACC@0.9', 'ACC_0.9']
# y_axis = ["SR2", "R2"]

# Define the data
data = {
    "Model": ["Qwen3 0.6B", "Qwen3 0.6B", "Qwen3 1.7B", "Qwen3 1.7B", 
              "Qwen3 4B", "Qwen3 4B", "Qwen3 8B", "Qwen3 8B", 
              "Qwen3 14B", "Qwen3 14B", "OpenAI o4-mini", "GPT-OSS-20B", "GPT-OSS-120B", "DeepSeek-V3-671B", "DeepSeek-R1-70B"],
    "Thinking": ["wo", "w", "wo", "w", "wo", "w", "wo", "w", "wo", "w", "w", "w", "w", "w", "w"],
    "R2": [0.739099, 0.749758, 0.474792, 0.799915, 
           0.657241, 0.713099, 0.761490, 0.799007, 
           0.611684, 0.591424, 0.865433, 0.803057, 0.878515, 0.846175, 0.874713772],
    "ACC_0.9": [0.692308, 0.730769, 0.346154, 0.769231,
                0.653846, 0.653846, 0.769231, 0.807692, 
                0.576923, 0.461538, 0.884615, 0.807692, 0.884615, 0.846154, 0.862068966],
    "total_tokens": [387606, 4642915, 621805, 1877705, 
                     543387, 7948242, 564913, 6106524, 
                     466224, 7039153, 2872857, 30000000, 24000000, 2400000, 3400000] 
}

df = pd.DataFrame(data)

# Add log of total tokens
df["log_tokens"] = np.log10(df["total_tokens"])

# Define marker for each model
model_markers = {
    "Qwen3 0.6B": "o",
    "Qwen3 1.7B": "s",
    "Qwen3 4B": "D",
    "Qwen3 8B": "^",
    "Qwen3 14B": "v",
    "OpenAI o4-mini": "P",
    "GPT-OSS-20B": "h",
    "GPT-OSS-120B": "H",
    "DeepSeek-V3-671B": "*",
    "DeepSeek-R1-70B": "X"
}

# Define colors for thinking
thinking_colors = {
    "w": "blue",    # with thinking
    "wo": "red"     # without thinking
}

# Define model sizes (in parameters) for scatter point sizing
model_sizes = {
    "Qwen3 0.6B": 0.6,
    "Qwen3 1.7B": 1.7,
    "Qwen3 4B": 4.0,
    "Qwen3 8B": 8.0,
    "Qwen3 14B": 14.0,
    "OpenAI o4-mini": 8.0,  # Estimated size
    "GPT-OSS-20B": 20.0,
    "GPT-OSS-120B": 120.0,
    "DeepSeek-V3-671B": 671.0,
    "DeepSeek-R1-70B": 70.0
}

# Function to convert model size to scatter point size
def get_scatter_size(model_size):
    # Scale the model size to appropriate scatter point sizes
    # Using larger base size and bigger multiplier for more dramatic differences
    base_size = 100
    size_multiplier = model_size *5
    return base_size + size_multiplier

# Linear regression on all points
X = df["log_tokens"].values.reshape(-1, 1)
y = df[y_axis[1]].values
reg = LinearRegression()
reg.fit(X, y)
x_range = np.linspace(df["log_tokens"].min(), df["log_tokens"].max(), 100).reshape(-1, 1)
y_pred = reg.predict(x_range)

# Plotting
fig, ax = plt.subplots(figsize=(9, 6))

# Plot each point with shape, color, and size based on model size
for _, row in df.iterrows():
    model = row["Model"]
    thinking = row["Thinking"]
    model_size = model_sizes.get(model, 1.0)  # Default to 1.0 if model not found
    scatter_size = get_scatter_size(model_size)
    
    ax.scatter(row["log_tokens"], row[y_axis[1]],
               marker=model_markers.get(model, "o"),  # Default to circle if marker not found
               color=thinking_colors[thinking],
               s=scatter_size,
               alpha=0.7)

# Regression line
ax.plot(x_range, y_pred, color='green', linestyle='--', label='Linear Fit')

# Custom legend
legend_elements = [
    # Individual models with their specific markers
    Line2D([0], [0], marker='o', color='w', label='Qwen3 0.6B', markerfacecolor='gray', markersize=10, alpha=0.7),
    Line2D([0], [0], marker='s', color='w', label='Qwen3 1.7B', markerfacecolor='gray', markersize=10, alpha=0.7),
    Line2D([0], [0], marker='D', color='w', label='Qwen3 4B', markerfacecolor='gray', markersize=10, alpha=0.7),
    Line2D([0], [0], marker='^', color='w', label='Qwen3 8B', markerfacecolor='gray', markersize=10, alpha=0.7),
    Line2D([0], [0], marker='v', color='w', label='Qwen3 14B', markerfacecolor='gray', markersize=10, alpha=0.7),
    Line2D([0], [0], marker='P', color='w', label='OpenAI o4-mini', markerfacecolor='gray', markersize=10, alpha=0.7),
    Line2D([0], [0], marker='h', color='w', label='GPT-OSS-20B', markerfacecolor='gray', markersize=10, alpha=0.7),
    Line2D([0], [0], marker='H', color='w', label='GPT-OSS-120B', markerfacecolor='gray', markersize=10, alpha=0.7),
    Line2D([0], [0], marker='*', color='w', label='DeepSeek-V3-671B', markerfacecolor='gray', markersize=10, alpha=0.7),
    Line2D([0], [0], marker='X', color='w', label='DeepSeek-R1-70B', markerfacecolor='gray', markersize=10, alpha=0.7),
    # Thinking colors
    Line2D([0], [0], marker='o', color='w', label='with thinking', markerfacecolor='blue', markersize=10, alpha=0.7),
    Line2D([0], [0], marker='o', color='w', label='without thinking', markerfacecolor='red', markersize=10, alpha=0.7),
    # Regression line
    Line2D([0], [0], color='green', linestyle='--', label='Linear Fit')
]
# bbox_to_anchor=(1.05, 1), loc='upper left'
ax.legend(loc='best', handles=legend_elements, fontsize=14)
ax.set_xlabel("log10(Generated Tokens)", fontsize=18)
ax.set_ylabel(y_axis[0], fontsize=18)
# ax.set_title("DE vs. log(Total Tokens) with Linear Fit")
ax.grid(True)

ax.tick_params(axis='both', which='major', labelsize=14)
ax.tick_params(axis='both', which='minor', labelsize=14)

plt.tight_layout()
plt.savefig("summary/compute_DE.pdf", dpi=300)
plt.savefig("summary/compute_DE.png", dpi=300)
plt.show()