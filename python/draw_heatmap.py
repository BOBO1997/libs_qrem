import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def draw_heatmap(data, row_labels=None, column_labels=None, norm=Normalize(vmin=-1, vmax=1)):

    if row_labels is None:
        row_labels = list(range(len(data)))
    if column_labels is None:
        column_labels = list(range(len(data[0])))
    # drawing part
    fig, ax = plt.subplots()
    # heatmap = ax.pcolor(data, cmap=plt.cm.Blues)
    heatmap = ax.pcolor(data, cmap="bwr", norm=norm)
    fig.colorbar(heatmap, ax=ax)

    ax.set_xticks(np.arange(data.shape[0]) + 0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1]) + 0.5, minor=False)

    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(row_labels, minor=False)
    ax.set_yticklabels(column_labels, minor=False)
    plt.show()
    # plt.savefig('image.png')

    return heatmap
