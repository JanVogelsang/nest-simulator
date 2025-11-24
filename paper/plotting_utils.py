import matplotlib
import matplotlib.pyplot as plt

single_column_in = 3.348
two_column_in = 7.09

matplotlib.use("pgf")
matplotlib.rcParams.update(
    {
        "pgf.texsystem": "pdflatex",
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
    }
)


def set_font_sizes(small=8, medium=10, large=12, family="Arial"):
    # plt.rc('text', usetex=True)
    plt.rc("font", size=small)  # controls default text sizes
    plt.rc("axes", titlesize=small)  # fontsize of the axes title
    plt.rc("axes", labelsize=medium)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=small)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=small)  # fontsize of the tick labels
    plt.rc("legend", fontsize=small)  # legend fontsize
    plt.rc("figure", titlesize=large)  # fontsize of the figure title
    plt.rc("font", family=family)


def save_grayscale(filename):
    from PIL import Image

    img = Image.open(filename).convert("L")
    name, ending = filename.split(".")
    img.save(name + "_g." + ending)
