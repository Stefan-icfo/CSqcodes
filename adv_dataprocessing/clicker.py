import matplotlib.pyplot as plt
import numpy as np

# Example raw data and axes (replace these with your actual data)
raw_data = np.random.rand(100, 100)  # This is a placeholder for your 2D array data
x_axis = np.linspace(-1.95, -1.85, raw_data.shape[1])  # Placeholder for your x-axis data
y_axis = np.linspace(-1.95, -1.85, raw_data.shape[0])  # Placeholder for your y-axis data

# Create a figure and display the raw data
fig, ax = plt.subplots()
# Use pcolormesh or imshow depending on how you want to visualize the raw data
c = ax.pcolormesh(x_axis, y_axis, raw_data, cmap='plasma')

# Add a colorbar to show the scale of the values
plt.colorbar(c, ax=ax)

# Lists to store the clicked points
clicks = []

# Define a function to capture clicks
def onclick(event):
    # Get the x and y pixel coordinates of the click
    x, y = event.xdata, event.ydata
    if x is not None and y is not None:
        clicks.append((x, y))
        print(f"Clicked on: ({x:.2f}, {y:.2f})")
        # Mark the clicked points
        ax.plot(x, y, 'ro')  # Red dot for marking the clicked point
        fig.canvas.draw()

# Connect the click event to the figure
cid = fig.canvas.mpl_connect('button_press_event', onclick)

# Display the image and wait for clicks
plt.show()

# After clicking, print the selected points
print("Selected points:", clicks)
