"""
This script generates barplots visualizing population distribution data from a 
PopulationDistribution.txt file.

It creates a barplot for each guide, dividing targets by the number of bulges 
and mismatches, and saves the resulting plots as PNG files. The script is intended 
to be run as a standalone program and expects four command-line arguments: the 
population distribution file, the maximum number of mismatches/bulges, the guide 
sequence, and a score or identifier for the output file.

Example of PopulationDistribution file (with 1 bulge max):
-Summary_CCATCGGTGGCCGTTTGCCCNNN
EAS     0,0     0,0     0,0     0,0     0,4     0,22    0,0     0,0     0,0     0,0
EUR     0,0     0,0     0,0     0,0     2,1     0,28    0,0     0,0     0,0     0,0
AFR     0,0     1,0     0,1     0,0     1,6     0,42    0,0     0,0     0,0     0,0
AMR     0,0     0,0     0,0     0,0     2,4     0,30    0,0     0,0     0,0     0,0
SAS     0,0     0,0     0,0     0,0     1,1     0,34    0,0     0,0     0,0     0,0
-Summary_GAGTCCGAGCAGAAGAAGAANNN
EAS     0,0     0,0     0,0     1,1     1,18    0,147   0,0     0,0     0,0     0,0
EUR     0,0     0,0     0,0     1,1     3,20    0,128   0,0     0,0     0,0     0,0
AFR     0,0     0,0     0,0     1,1     3,31    0,253   0,0     0,0     0,0     0,0
AMR     0,0     0,0     0,0     1,1     3,26    0,146   0,0     0,0     0,0     0,0
SAS     0,0     0,0     0,0     2,1     3,21    0,153   0,0     0,0     0,0     0,0
"""

from typing import Tuple, Optional, Dict, List
from matplotlib.legend_handler import HandlerTuple


import matplotlib.pyplot as plt
import matplotlib.colors as mc
import numpy as np

import matplotlib
import colorsys
import warnings
import math
import sys


warnings.filterwarnings("ignore")  # suppress warnings

TARGETSCOLORS = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c"]
TITLESIZE = 18
FONTSIZE = 17

# set matplotlib params
matplotlib.use("Agg")  # do not use X11
plt.rcParams["figure.dpi"] = 400  # set matplotlib for pdf editing
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.style.use("seaborn-poster")


def adjust_lightness(
    color: str, amount: Optional[float] = 0.5
) -> Tuple[float, float, float]:
    """Adjusts the lightness of a given color by a specified amount.

    This function converts the input color to HLS, modifies its lightness, and
    returns the adjusted RGB value.

    Args:
        color: The color to adjust, as a string (name or hex code).
        amount: The factor by which to adjust the lightness (default is 0.5).

    Returns:
        A tuple representing the adjusted RGB color.

    Raises:
        ValueError: If the color cannot be converted to RGB.
    """
    assert amount is not None
    try:
        c = mc.cnames[color]
    except KeyError:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def _process_data_line(
    line: str,
    mmbul_num: int,
    barplot_data: Dict[str, Dict[int, List[int]]],
    current_max: int,
    number_bars: int,
) -> Tuple[int, int]:
    """Parses a line of population distribution data and updates the barplot data.

    This function extracts population counts from a line, updates the barplot
    data structure, and tracks the maximum value and number of bars.

    Args:
        line: The input line containing population data.
        mmbul_num: The maximum number of mismatches or bulges.
        barplot_data: The dictionary to update with parsed data.
        current_max: The current maximum value found.
        number_bars: The current number of bars.

    Returns:
        A tuple containing the updated maximum value and number of bars.

    Raises:
        ValueError: If the line is malformed or contains invalid data.
    """
    if not line.strip():
        return current_max, number_bars
    fields = line.split()
    if len(fields) < mmbul_num + 2:  # population name + (mmbul_num + 1) data fields
        raise ValueError(
            f"Insufficient fields in line: expected {mmbul_num + 2}, got {len(fields)}"
        )
    population_name = fields[0]
    if not population_name or population_name.isspace():  # validate population name
        raise ValueError(
            f"Empty or whitespace-only population name? ({population_name})"
        )
    barplot_data[population_name] = {}  # initialize population data
    line_max_value = current_max
    for i in range(mmbul_num + 1):  # process each count field
        try:  # parse comma-separated values
            values = [
                int(x.strip()) for x in fields[i + 1].split(",") if x.strip()
            ] or [
                0
            ]  # default to 0 if no valid values
            number_bars = len(fields[i + 1].split(","))
            barplot_data[population_name][i] = values
            line_max_value = max(line_max_value, sum(values))
        except ValueError as e:
            raise ValueError(
                f"Invalid integer values in field '{fields[i + 1]}'"
            ) from e
    return line_max_value, number_bars


def read_population_distribution(
    popdist_fname: str, guide: str, mmbul_num: int
) -> Tuple[Dict[str, Dict[int, List[int]]], int, int]:
    """Reads a population distribution file and extracts barplot data for a
    specific guide.

    This function parses the population distribution file, finds the section for
    the given guide, and collects the relevant data for plotting. It returns the
    barplot values, the maximum value found, and the number of bars.

    Args:
        popdist_fname: Path to the population distribution file.
        guide: The guide sequence to search for in the file.
        mmbul_num: The maximum number of mismatches or bulges.

    Returns:
        A tuple containing the barplot values dictionary, the maximum value, and
            the number of bars.

    Raises:
        FileNotFoundError: If the population distribution file does not exist.
        IOError: If there is an error reading the file.
        ValueError: If the file contains malformed lines or inconsistent data.
    """
    barplot_values = {}
    max_value = 0
    guide_found = False
    number_bars = -1
    number_bars_computed = False
    try:
        with open(popdist_fname, "r", encoding="utf-8") as infile:
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                if guide in line:
                    guide_found = True
                    continue
                if not guide_found:
                    continue
                if "-Summary_" in line:
                    break
                try:
                    max_value, number_bars_ = _process_data_line(
                        line, mmbul_num, barplot_values, max_value, number_bars
                    )
                    if not number_bars_computed and number_bars_ != number_bars:
                        number_bars = number_bars_
                    if number_bars_computed and number_bars_ != number_bars:
                        raise ValueError(
                            f"Mismatching number of bars at line {line}, expected {number_bars}, got {number_bars_}"
                        )
                except (ValueError, IndexError) as e:
                    print(f"Warning: Line {line_num}: {e}")
                    continue
    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found: {popdist_fname}") from e
    except IOError as e:
        raise IOError(f"Error reading file {popdist_fname}: {e}") from e
    return barplot_values, max_value, number_bars


def _calculate_population_sum(
    pop_data: Dict[int, List[int]], mmbul_num: int, number_bars: int
) -> int:
    """Calculates the sum of population counts across all categories up to a
    given threshold.

    This function aggregates the counts for each bar position across all mismatch
    or bulge categories below the maximum.

    Args:
        pop_data: Dictionary mapping category indices to lists of counts.
        mmbul_num: The maximum number of mismatches or bulges to include.
        number_bars: The number of bar positions to sum over.

    Returns:
        The total sum of counts across all specified categories and bar positions.
    """
    bar_sums = [0] * number_bars  # initialize sum array for each bar position
    # sum across all count categories (0 to mmbul_num-1)
    for i in range(mmbul_num):
        for bar_idx, value in enumerate(pop_data[i]):
            if bar_idx < number_bars:  # safety check
                bar_sums[bar_idx] += value
    return sum(bar_sums)


def compute_lower_bar_values(
    barplot_values: Dict[str, Dict[int, List[int]]], mmbul_num: int, number_bars: int
) -> Dict[str, int]:
    """Computes the sum of lower bar values for each population.

    This function calculates the total count for each population by summing all
    values from count indices below the maximum threshold.

    Args:
        barplot_values: Dictionary mapping population names to their count data.
        mmbul_num: The maximum number of mismatches or bulges to include.
        number_bars: The number of bar positions to sum over.

    Returns:
        A dictionary mapping each population to its total lower bar value.
    """
    # calculate lower bound sums for each population
    lower_barplot_values = {}
    for population, pop_data in barplot_values.items():
        # Sum all values from count indices 0 to total-1
        population_total = _calculate_population_sum(pop_data, mmbul_num, number_bars)
        lower_barplot_values[population] = population_total
    return lower_barplot_values


def _create_default_range() -> Tuple[np.ndarray, bool]:
    """Creates a default y-axis range for plots when no data is available.

    This function returns a default numpy array and a flag indicating no data
    was found.

    Returns:
        A tuple containing a default numpy array for the y-axis and a boolean
            flag set to True.
    """
    return np.arange(0, 1, 1), True


def _create_dynamic_range(total_max: int) -> np.ndarray:
    """Generates a dynamic y-axis range based on the maximum value for plotting.

    This function calculates a step size and upper limit to create a suitable
    numpy array for the y-axis.

    Args:
        total_max: The maximum value to base the y-axis range on.

    Returns:
        A numpy array representing the dynamic y-axis range.
    """
    # calculate step size (20% of max for better granularity)
    step_size = max(1, math.ceil(total_max / 5))
    # calculate upper limit with 10% padding
    upper_limit = total_max + math.ceil(total_max / 10)
    # create range starting from 0
    return np.arange(0, upper_limit + step_size, step_size)


def _compute_plot_yrange(
    lower_barplot_values: Dict[str, int], max_value: int
) -> Tuple[np.ndarray, bool]:
    """Computes the y-axis range for a barplot based on lower bar values and the
    maximum value.

    This function determines an appropriate y-axis range for plotting, using the
    sum of lower bar values and the provided maximum.

    Args:
        lower_barplot_values: Dictionary mapping populations to their lower bar
            values.
        max_value: The maximum value found in the data.

    Returns:
        A tuple containing a numpy array for the y-axis range and a boolean
            indicating if no data was found.
    """
    if not lower_barplot_values:  # validate inputs
        return _create_default_range()
    # find maximum value across all lower bounds
    max_lower_value = max(lower_barplot_values.values())
    total_max = max_value + max_lower_value  # update total maximum
    if total_max <= 0:  # calculate appropriate Y-range
        return _create_default_range()
    y_range = _create_dynamic_range(total_max)
    return y_range, False


def compute_plot_yrange(
    lower_barplot_values: Dict[str, int], max_value: int
) -> Tuple[np.ndarray, bool]:
    """Safely computes the y-axis range for a barplot using lower bar values and
    the maximum value.

    This function wraps the internal y-range computation and handles errors by
    returning a default range if needed.

    Args:
        lower_barplot_values: Dictionary mapping populations to their lower bar
            values.
        max_value: The maximum value found in the data.

    Returns:
        A tuple containing a numpy array for the y-axis range and a boolean
            indicating if no data was found.
    """
    try:
        return _compute_plot_yrange(lower_barplot_values, max_value)
    except (ValueError, TypeError, KeyError) as e:
        # Handle specific expected errors
        print(f"Warning: Error calculating Y-range: {e}")
        return _create_default_range()
    except Exception as e:
        # Handle unexpected errors
        print(f"Unexpected error in Y-range calculation: {e}")
        return _create_default_range()


def assign_targets_colors(number_bars: int) -> List[str]:
    """Assigns a list of colors for barplot targets based on the number of bars.

    This function returns a list of color codes, extending or trimming the
    default color set as needed to match the number of bars.

    Args:
        number_bars: The number of bars to assign colors to.

    Returns:
        A list of color codes for the barplot.
    """
    colors = TARGETSCOLORS.copy()
    if len(colors) < number_bars:  # must extend target colors
        for _, c in mc.TABLEAU_COLORS.items():  # type: ignore
            colors.append(c)
            if len(colors) == number_bars:
                break
    elif number_bars < len(colors):  # can trim target colors
        colors = colors[:number_bars]
    return colors


def retrieve_lower_counts(lower_barplot_values: Dict[str, int]) -> List[int]:
    """Retrieves the lower bar values as a list of counts for each population.

    This function extracts the total lower bar value for each population from
    the input dictionary.

    Args:
        lower_barplot_values: Dictionary mapping populations to their lower bar
            values.

    Returns:
        A list of lower bar values for all populations.
    """
    return [v for _, v in lower_barplot_values.items()]


def retrieve_bar_counts(
    barplot_values: Dict[str, List[int]], number_bars: int
) -> List[List[int]]:
    """Retrieves the counts for each bar position across all populations.

    This function returns a list of lists, where each inner list contains the
    counts for a specific bar position for all populations.

    Args:
        barplot_values: Dictionary mapping populations to their bar counts.
        number_bars: The number of bar positions to retrieve.

    Returns:
        A list of lists, each containing the counts for a bar position across
            all populations.
    """
    return [[v[i] for v in barplot_values.values()] for i in range(number_bars)]


def retrieve_barplot_counts(
    barplot_values: Dict[str, List[int]],
    lower_barplot_values: Dict[str, int],
    number_bars: int,
) -> Tuple[List[int], List[List[int]]]:
    """Retrieves both lower bar values and per-bar counts for all populations.

    This function returns a tuple containing the list of lower bar values and a
    list of lists with counts for each bar position.

    Args:
        barplot_values: Dictionary mapping populations to their bar counts.
        lower_barplot_values: Dictionary mapping populations to their lower bar
            values.
        number_bars: The number of bar positions to retrieve.

    Returns:
        A tuple with a list of lower bar values and a list of lists for bar counts.
    """
    lower_counts = retrieve_lower_counts(lower_barplot_values)
    bar_counts = retrieve_bar_counts(barplot_values, number_bars)
    return lower_counts, bar_counts


def compute_legend_labels_n_colors(
    barplot_values: Dict[str, List[int]],
    bar_low: List,
    bars: List,
    mmbul_num: int,
    lower_vals: bool,
    number_bars: int,
) -> Tuple[List[str], List[int]]:
    """Generates legend labels and color handles for the barplot.

    This function creates descriptive labels and corresponding color handles for
    each bar in the plot, based on the mismatch and bulge counts.

    Args:
        barplot_values: Dictionary mapping populations to their bar counts.
        bar_low: List containing the lower bar plot objects.
        bars: List of bar plot objects for each bar position.
        mmbul_num: The maximum number of mismatches or bulges.
        lower_vals: Boolean indicating if lower bar values are present.
        number_bars: The number of bar positions.

    Returns:
        A tuple containing a list of legend labels and a list of color handles.
    """
    labels = [f"MM+B < {mmbul_num}"] if lower_vals else []
    handles_colors = [bar_low[0][0]] if lower_vals else []
    for i in range(number_bars):
        if i == 0 and barplot_values["REF"][0] != 0:
            labels.append(f"{mmbul_num} MM")  # mms + 0 bulges
            handles_colors.append(bars[0][0])
        elif mmbul_num - i < 0 and barplot_values["REF"][i] != 0:
            # avoid negative numbers on MM values
            labels.append(f"0 MM + {i} B")
            handles_colors.append(bars[i][0])
        elif mmbul_num - i >= 0 and barplot_values["REF"][i]:
            labels.append(f"{mmbul_num - i} MM + {i} B")
            handles_colors.append(bars[i][0])
    # reverse lists
    labels.reverse()
    handles_colors.reverse()
    return labels, handles_colors


def draw_barplot(
    labels: List[str],
    barplot_values: Dict[str, List[int]],
    lower_barplot_values: Dict[str, int],
    max_value: int,
    mmbul_num: int,
    number_bars: int,
    guide: str,
    score: str,
) -> None:
    """Draws and saves a barplot visualizing population distribution data.

    This function creates a barplot for the given guide and population data,
    displaying the distribution of targets by mismatches and bulges, and saves
    the plot as a PNG file.

    Args:
        labels: List of population labels for the x-axis.
        barplot_values: Dictionary mapping populations to their bar counts.
        lower_barplot_values: Dictionary mapping populations to their lower bar
            values.
        max_value: The maximum value found in the data.
        mmbul_num: The maximum number of mismatches or bulges.
        number_bars: The number of bar positions.
        guide: The guide sequence used for the plot.
        score: The score or identifier for the plot file name.

    Returns:
        None. The function saves the generated barplot as a PNG file.
    """
    ax = plt.figure()  # initialize plot
    # compute plot x-range and y-range
    x_range = np.arange(0, len(labels), 1)
    y_range, nodata = compute_plot_yrange(lower_barplot_values, max_value)
    colors = assign_targets_colors(number_bars)  # retrieve barplot colors
    width = 0.5  # bars width
    lower_counts, bar_counts = retrieve_barplot_counts(
        barplot_values, lower_barplot_values, number_bars
    )  # retrieve bar counts
    bar_low = []  # bars storing the counts of targets lower than threshold
    if max(lower_counts) > 0:  # lower targets found
        bar_low = [
            plt.bar(
                x_range,
                lower_counts,
                width=width,
                color="yellow",
                align="edge",
                edgecolor="black",
            )
        ]
    bars, previous_bar = [], lower_counts
    for i in range(number_bars):  # format counts for targets above threshold
        bars.append(
            plt.bar(
                x_range,
                bar_counts[i],
                width=width,
                color=colors[i],
                align="edge",
                bottom=previous_bar,
                edgecolor="black",
            )
        )
        previous_bar = [x + previous_bar[j] for j, x in enumerate(bar_counts[i])]
    legend_labs, handles_colors = compute_legend_labels_n_colors(
        barplot_values, bar_low, bars, mmbul_num, max(lower_counts) > 0, number_bars
    )  # compute legened labels and color handles, and create legend
    plt.legend(handles_colors, legend_labs, fontsize=FONTSIZE, handlelength=5, handler_map={tuple: HandlerTuple(ndivide=None)}, title="MM mismatches, B bulges", title_fontsize=FONTSIZE)  # type: ignore
    plt.title(
        f"Targets with up to {mmbul_num} mismatches and/or bulges by "
        "superpopulation",
        size=TITLESIZE,
    )
    if nodata:  # no barplot data to display -> show message
        plt.annotate(
            f"No targets found with {mmbul_num} mismatches + bulges",
            (2.5, 0),
            size=FONTSIZE,
            ha="center",
            va="center",
        )
        return
    plt.xticks(x_range + 0.25, labels, size=FONTSIZE)  # style x axis
    digits = int(math.log10(max_value)) + 1
    size_y_ticks = max(16, FONTSIZE - (2 * (digits - 5))) if digits > 5 else FONTSIZE
    plt.yticks(y_range, size=size_y_ticks)  # style y axis
    plt.tight_layout()
    plt.savefig(
        f"populations_distribution_{guide}_{mmbul_num}total_{score}_new.png",
        format="png",
    )


def create_population_dist_plot() -> None:
    """Creates and saves a population distribution barplot from a
    PopulationDistribution.txt file.

    This function parses command-line arguments, reads the population distribution
    data for a specific guide, computes the necessary values, and generates a
    barplot image.

    Returns:
        None. The function saves the generated barplot as a PNG file.
    """
    assert len(sys.argv[1:]) == 4
    popdist_fname, mmbul_num, guide, score = sys.argv[1:]  # read input args
    mmbul_num = int(mmbul_num)
    # parse population distribution file for input guide
    barplot_values, max_value, number_bars = read_population_distribution(
        popdist_fname, guide, mmbul_num
    )
    # compute lower bars values
    lower_barplot_values = compute_lower_bar_values(
        barplot_values, mmbul_num, number_bars
    )
    # extract final values (at index 'mmbul_num')
    final_barplot_values = {
        pop: data[mmbul_num] for pop, data in barplot_values.items()
    }
    # draw barplot
    populations = list(final_barplot_values.keys())
    draw_barplot(
        populations,
        final_barplot_values,
        lower_barplot_values,
        max_value,
        mmbul_num,
        number_bars,
        guide,
        score,
    )


if __name__ == "__main__":
    create_population_dist_plot()