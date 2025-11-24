from collections import namedtuple

import pandas as pd
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
from plotting_utils import *

size_factor = 1.3
matplotlib_params = {
    "text.latex.preamble": ["\\usepackage{gensymb}"],
    "image.origin": "lower",
    "image.interpolation": "nearest",
    "axes.grid": False,
    "axes.labelsize": 15 * size_factor,
    "axes.titlesize": 19 * size_factor,
    "font.size": 16 * size_factor,
    "legend.fontsize": 11 * size_factor,
    "xtick.labelsize": 11 * size_factor,
    "ytick.labelsize": 11 * size_factor,
    "text.usetex": False,
}

additional_params = {"figsize_single": [6.1 * 1.5, 6.1], "figsize_double": [12.2, 6.1 * 1.1]}

light = namedtuple(
    "light",
    [
        "light_blue",
        "orange",
        "light_yellow",
        "pink",
        "light_cyan",
        "mint",
        "pear",
        "olive",
        "pale_grey",
        "black",
    ],
)(
    light_blue="#77AADD",
    orange="#EE8866",
    light_yellow="#EEDD88",
    pink="#FFAABB",
    light_cyan="#99DDFF",
    mint="#44BB99",
    pear="#BBCC33",
    olive="#AAAA00",
    pale_grey="#DDDDDD",
    black="#000000",
)

color_params = {
    "wall_time_total": light.pale_grey,
    "sim_factor": light.pink,
    "phase_total_factor": light.orange,
    "time_simulate": light.pink,
    "time_construction_create+time_construction_connect": light.light_cyan,
    "time_update_spike_data": light.orange,
    "time_deliver_spike_data": light.light_blue,
    "time_communicate_spike_data": light.mint,
    "wall_time_phase_collocate": light.light_yellow,
    "frac_phase_update": light.orange,
    "frac_phase_deliver": light.light_blue,
    "frac_phase_communicate": light.mint,
    "frac_phase_collocate": light.light_yellow,
    "phase_update_factor": light.orange,
    "phase_deliver_factor": light.light_blue,
    "phase_communicate_factor": light.mint,
    "phase_collocate_factor": light.light_yellow,
    "total_memory": light.olive,
    "total_memory_per_node": light.pear,
    "phase_omp_sync_factor": light.olive,
    "phase_mpi_sync_factor": light.pink,
    "frac_phase_omp_sync": light.olive,
    "frac_phase_mpi_sync": light.pink,
}

label_params = {
    "threads_per_node": "OMP threads per node",
    "tasks_per_node": "MPI processes per node",
    "num_nodes": "Nodes",
    "num_ranks": "MPI processes",
    "wall_time_total": "Total",
    "wall_time_preparation": "Preparation",
    "wall_time_presim": "Presimulation",
    "wall_time_creation": "Creation",
    "wall_time_connect": "Connection",
    "time_simulate": "State propagation",
    "wall_time_phase_total": "All phases",
    "time_update_spike_data": "Update",
    "wall_time_phase_collocate": "Collocation",
    "time_communicate_spike_data": "Communication",
    "time_deliver_spike_data": "Delivery",
    "time_construction_create+time_construction_connect": "Network construction",
    "max_memory": "Memory",
    "sim_factor": "State propagation",
    "frac_phase_update": "Update",
    "frac_phase_communicate": "Communication",
    "frac_phase_deliver": "Delivery",
    "frac_phase_collocate": "Collocation",
    "phase_update_factor": "Update",
    "phase_communicate_factor": "Communication",
    "phase_deliver_factor": "Delivery",
    "phase_collocate_factor": "Collocation",
    "phase_total_factor": "All phases",
    "total_memory": "Memory",
    "total_memory_per_node": "Memory per node",
    "phase_omp_sync_factor": "Thread synchronization",
    "phase_mpi_sync_factor": "Mpi synchronization",
    "frac_phase_omp_sync": "Thread synchronization",
    "frac_phase_mpi_sync": "Mpi synchronization",
}


class Plot:
    def __init__(
        self,
        x_axis,
        x_ticks="data",
        data_file="/path/to/data",
        time_scaling=1,
        default_use_inter_area_axonal_delay=False,
    ):
        self.x_axis = x_axis
        self.x_ticks = x_ticks
        self.matplotlib_params = matplotlib_params
        self.additional_params = additional_params
        self.color_params = color_params
        self.label_params = label_params
        self.time_scaling = time_scaling

        self.load_data(data_file, default_use_inter_area_axonal_delay)
        self.compute_derived_quantities()

    def load_data(self, data_file, default_use_inter_area_axonal_delay):
        self.df = pd.read_csv(data_file, delimiter=",")

        if default_use_inter_area_axonal_delay:
            self.df["use_inter_area_axonal_delay"] = False

        self.df["num_ranks"] = self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["time_update_spike_data"] = self.df["time_update_mean"]
        self.df["time_deliver_spike_data"] = self.df["time_deliver_spike_data_mean"]
        self.df["time_omp_sync_spike_data"] = self.df["time_omp_synchronization_simulation"]
        self.df["time_mpi_sync_spike_data"] = self.df["time_mpi_synchronization"]
        self.df["time_communicate_spike_data"] = (
            self.df["time_communicate_spike_data"] - self.df["time_mpi_synchronization"]
        )

        self.df["base_memory"] /= self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["network_memory"] /= self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["init_memory"] /= self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["total_memory"] /= self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["base_memory_rss"] /= self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["network_memory_rss"] /= self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["init_memory_rss"] /= self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["total_memory_rss"] /= self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["base_memory_peak"] /= self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["network_memory_peak"] /= self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["init_memory_peak"] /= self.df["num_nodes"] * self.df["tasks_per_node"]
        self.df["total_memory_peak"] /= self.df["num_nodes"] * self.df["tasks_per_node"]

        dict_ = {
            "num_nodes": "first",
            "threads_per_task": "first",
            "tasks_per_node": "first",
            "num_ranks": "first",
            "model_time_sim": "first",
            "network_state": "first",
            "use_inter_area_axonal_delay": "first",
            "time_construction_create": ["mean", "std"],
            "time_construction_connect": ["mean", "std"],
            "time_simulate": ["mean", "std"],
            "time_communicate_prepare": ["mean", "std"],
            "py_time_create": ["mean", "std"],
            "py_time_connect": ["mean", "std"],
            "py_time_presimulate": ["mean", "std"],
            "py_time_simulate": ["mean", "std"],
            "base_memory": ["mean", "std"],
            "network_memory": ["mean", "std"],
            "init_memory": ["mean", "std"],
            "total_memory": ["mean", "std"],
            "base_memory_rss": ["mean", "std"],
            "network_memory_rss": ["mean", "std"],
            "init_memory_rss": ["mean", "std"],
            "total_memory_rss": ["mean", "std"],
            "base_memory_peak": ["mean", "std"],
            "network_memory_peak": ["mean", "std"],
            "init_memory_peak": ["mean", "std"],
            "total_memory_peak": ["mean", "std"],
            "num_connections": ["mean", "std"],
            "local_spike_counter": ["mean", "std"],
            "time_collocate_spike_data": ["mean", "std"],
            "time_communicate_spike_data": ["mean", "std"],
            "time_deliver_spike_data": ["mean", "std"],
            "time_communicate_target_data": ["mean", "std"],
            "time_gather_spike_data": ["mean", "std"],
            "time_gather_target_data": ["mean", "std"],
            "time_omp_synchronization_construction": ["mean", "std"],
            "time_omp_synchronization_simulation": ["mean", "std"],
            "time_update_spike_data": ["mean", "std"],
            "time_mpi_sync_spike_data": ["mean", "std"],
            "time_omp_sync_spike_data": ["mean", "std"],
        }

        self.df = (
            self.df.drop("rng_seed", axis=1)
            .groupby(
                [
                    "num_nodes",
                    "threads_per_task",
                    "tasks_per_node",
                    "model_time_sim",
                    "model_time_presim",
                    "use_inter_area_axonal_delay",
                    "network_state",
                ],
                as_index=False,
            )
            .agg(dict_)
        )

    def compute_derived_quantities(self):
        """
        Do computations to get parameters needed for plotting.
        """
        self.df["num_nvp"] = self.df["threads_per_task"] * self.df["tasks_per_node"]
        self.df["model_time_sim"] /= self.time_scaling
        self.df["sim_factor"] = self.df["time_simulate"].iloc[:, 0] / self.df["model_time_sim"].iloc[:, 0]

        self.df["time_phase_total"] = (
            self.df["time_update_spike_data"].iloc[:, 0]
            + self.df["time_communicate_spike_data"].iloc[:, 0]
            + self.df["time_deliver_spike_data"].iloc[:, 0]
            + self.df["time_collocate_spike_data"].iloc[:, 0]
        )
        self.df["phase_total_factor"] = self.df["time_phase_total"] / self.df["model_time_sim"].iloc[:, 0]

        for phase in ["update", "communicate", "deliver", "collocate", "omp_sync", "mpi_sync"]:
            self.df["phase_" + phase + "_factor"] = (
                self.df["time_" + phase + "_spike_data"].iloc[:, 0] / self.df["model_time_sim"].iloc[:, 0]
            )

            self.df["frac_phase_" + phase] = (
                100 * self.df["time_" + phase + "_spike_data"].iloc[:, 0] / self.df["time_phase_total"]
            )

        self.df["total_memory_per_node"] = self.df["total_memory"].iloc[:, 0] / self.df["num_nodes"].iloc[:, 0]


x_axis = "num_nodes"
y_axis = "sim_factor"
save_path = "MAM_Jureca_ax-delay.png"

B_correction = Plot(data_file="data/MAM_Jureca_correction.csv", x_axis=x_axis, time_scaling=int(1e3))
B_master = Plot(
    data_file="data/MAM_Jureca_master.csv",
    x_axis=x_axis,
    time_scaling=int(1e3),
    default_use_inter_area_axonal_delay=True,
)
B_correction.df = B_correction.df[B_correction.df["network_state"].iloc[:, 0] == "ground"]
B_master.df = B_master.df[B_master.df["network_state"].iloc[:, 0] == "ground"]


def plot_correction_vs_no_crrection():
    cmap = plt.get_cmap("tab10")
    fig, ax = plt.subplots(figsize=(single_column_in, 3.5))
    set_font_sizes()
    plt.setp(ax.spines.values(), linewidth=2)

    B = B_correction
    for i, label in enumerate(["With axonal delay", "Without axonal delay"]):
        colors = cm.get_cmap("Blues", 5) if i == 0 else cm.get_cmap("Reds", 5)
        df = B.df.loc[B.df["use_inter_area_axonal_delay"].iloc[:, 0] == (i == 0)]
        x_data = df[x_axis].iloc[:, 0]
        y_data = df[y_axis]
        ax.plot(x_data, y_data, label=label, color=colors(2), marker=("x" if i == 0 else "o"))

    ax.set_xlabel("Number of nodes")
    ax.set_ylabel("Simulation time [s]")
    ax.set_xscale("log")
    # ax.set_yscale("log")
    ax.minorticks_off()
    ax.tick_params(width=2)
    fig.tight_layout()
    fig.savefig("benchmark_axonal_vs_dendritic.pgf")
    fig.savefig("benchmark_axonal_vs_dendritic.jpg", dpi=300)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xticks(B.df[x_axis].iloc[:, 0])
    ax.set_xticklabels(B.df[x_axis].iloc[:, 0])
    plt.title("MAM, Jureca, axonal vs dendritic (inter-area)")
    plt.legend()
    plt.show()
    plt.close()


def plot_ax_vs_master():
    fig, ax = plt.subplots(figsize=(single_column_in, 3.5))
    set_font_sizes()
    plt.setp(ax.spines.values(), linewidth=2)

    df_corr = B_correction.df.loc[B_correction.df["use_inter_area_axonal_delay"].iloc[:, 0] == False]
    df_master = B_master.df
    for i, [label, df] in enumerate([("Correction-based", df_corr), ("Master", df_master)]):
        colors = cm.get_cmap("Reds", 5) if i == 0 else cm.get_cmap("Oranges", 5)
        x_data = df[x_axis].iloc[:, 0]
        y_data = df[y_axis]
        ax.plot(x_data, y_data, label=label, color=colors(2), marker=("x" if i == 0 else "o"))

    ax.set_xlabel("Number of nodes")
    ax.set_ylabel("Real-time factor")
    ax.set_xscale("log")
    # ax.set_yscale("log")
    ax.minorticks_off()
    ax.tick_params(width=2)
    fig.tight_layout()
    fig.savefig("benchmark_correction_vs_master.pgf")
    fig.savefig("benchmark_correction_vs_master.jpg", dpi=300)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xticks(df_corr[x_axis].iloc[:, 0])
    ax.set_xticklabels(df_corr[x_axis].iloc[:, 0])
    plt.title("MAM, Jureca, fully dendritic, Correction vs Master")
    plt.legend()
    plt.show()
    plt.close()


plot_correction_vs_no_crrection()
plot_ax_vs_master()
