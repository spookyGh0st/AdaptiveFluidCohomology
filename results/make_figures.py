import fnmatch
import itertools
import math
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import cv2
import numpy as np
import glob
from concurrent.futures import ProcessPoolExecutor
from PIL import Image
# from scipy.signal import savgol_filter

mpl.use("Agg")

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 12,
})

LINE_STYLES = ["-", "--", "-.", ":"]
MARKERS = ["o", "s", "^", "d"]

STYLE_COMBINATIONS = [
    (ls, mk)
    for mk in MARKERS
    for ls in LINE_STYLES
]

def assign_colors(names):
    """
    Assign a unique color to each short_name across all plots.
    
    Args:
        results (list): list of short_name
    Returns:
        dict: mapping {short_name: color}
    """
    cmap = cm.get_cmap("Dark2", None)  # use tab10 categorical colormap
    return {name: cmap(i) for i, name in enumerate(names)}

def group_columns(columns):
    """
    Group columns by everything except IDX(n) and replace IDX(n) with _n for DataFrame access.
    """
    grouped = defaultdict(list)
    
    for col in columns:
        if col.lower() == "time":
            continue
        
        # Find all IDX(...) occurrences
        idx_matches = re.findall(r'IDX\((\d+)\)', col)
        if idx_matches:
            # Convert to integers for sorting
            idxs = [int(x) for x in idx_matches]
            # Group key: string with all IDX(...) removed
            key = re.sub(r'_?IDX\(\d+\)', '', col)
            # Replace IDX(n) with _{n} for DataFrame access
            new_col = re.sub(r'IDX\((\d+)\)', r'{\1}', col)
            grouped[key].append((col, new_col))
        else:
            # No IDX, group by full string
            grouped[col].append((col, col))
    
    # Sort items within each group by their numeric indices
    result = {}
    for key, items in grouped.items():
        sorted_items = [(coln,coll) for coln, coll in sorted(items, key=lambda x: x[0])]
        result[key] = sorted_items
    
    return result


def plot_grouped_columns(results, output_dir):
    if not results:
        print("No data to plot.")
        return
    
    # Determine all column groups from the first dataframe
    all_columns = results[0][2].columns
    grouped_columns = group_columns(all_columns)


    names = [name for name, _, _ in results]
    color_map = assign_colors(names)
    # Generate plots
    for group_name, cols in grouped_columns.items():
        plt.figure(figsize=(5, 3))

        for short_name, label, df in results:
            color = color_map[short_name]
            if len(cols) == 1:
                (coln,coll) = cols[0]
                plt.plot(df["time"], df[coln],
                         label=label,
                         color=color,
                         linestyle=LINE_STYLES[0])
            else:
                # Different line style for each column in group
                for style, (coln,coll) in zip(itertools.cycle(STYLE_COMBINATIONS), cols):
                    plt.plot(df["time"], df[coln],
                            label=f"{short_name} {coll}",
                            color=color,
                            linestyle=style[0],
                            marker = style[1],
                            )
        if len(cols) > 1: 
            custom_lines = []
            from matplotlib.lines import Line2D
            for style, (coln,coll) in zip(itertools.cycle(STYLE_COMBINATIONS), cols):
                custom_lines.append(
                    Line2D([0],[0],color="black",linestyle=style[0],label=coll,marker=style[1])
                )
            plt.legend(handles=custom_lines,fontsize="xx-small")
        # ax.set_xlabel("Time")
        # ax.set_ylabel(group_name)

        plt.xlabel("Time")
        plt.ylabel(group_name)
        plt.title(group_name)
        plt.legend()
        plt.tight_layout()

        outfile = os.path.join(output_dir, f"{group_name}.pdf")
        plt.savefig(outfile)
        plt.close()
        print(f"Saved {outfile}")

def plot_bigpage(results, output_file,config):
    """
    Create one combined PDF with all grouped plots.
    
    Args:
        results (list): (short_name, DataFrame) pairs
        output_file (str): path to save combined PDF
    """
    if not results:
        print("No data to plot.")
        return
    ncols, figsize = config

    # Determine groups
    all_columns = results[0][2].columns
    grouped_columns = group_columns(all_columns)

    # Prepare subplot grid
    n_groups = len(grouped_columns)
    nrows = math.ceil(n_groups / ncols)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=(figsize[0] * ncols, figsize[1] * nrows),
                             constrained_layout=True)

    # Flatten axes array for easy iteration
    if isinstance(axes, plt.Axes):
        axes = [axes]
    else:
        axes = axes.ravel()

    # Assign consistent colors
    names = [name for name, _, _ in results]
    color_map = assign_colors(names)

    for ax, (group_name, cols) in zip(axes, grouped_columns.items()):
        if len(cols) == 0 :
            continue
        for short_name,label, df in results:
            color = color_map[short_name]

            if len(cols) == 1:
                col = cols[0][0]
                # yhat = savgol_filter(df[col], 100, 3) 
                ax.plot(df["time"], df[col],
                        label=label,
                        color=color,
                        linestyle=LINE_STYLES[0])
            else:
                for (ls,ms), (coln,coll) in zip(itertools.cycle(STYLE_COMBINATIONS), cols):
                    # yhat = savgol_filter(df[coln], 100, 3) 
                    ax.plot(df["time"], df[coln],
                            label=f"{short_name} {coll}",
                            color=color,
                            linestyle=ls, 
                            markersize=2,
                            markevery=100,
                            marker=ms)
        if len(cols) > 1: 
            custom_lines = []
            from matplotlib.lines import Line2D
            for (ls,mk), (coln,coll) in zip(itertools.cycle(STYLE_COMBINATIONS), cols):
                custom_lines.append(
                    Line2D([0],[0],color="black",linestyle=ls, marker=mk,label=coll,
                            markersize=2,
                            markevery=25,
                           )
                )
            ax.legend(handles=custom_lines,fontsize="xx-small")
        # ax.set_xlabel("Time")
        # ax.set_ylabel(group_name)
        ax.set_title(group_name)

    # Add legend once (outside plots)
    handles, labels = axes[0].get_legend_handles_labels()
    # leg = fig.legend(handles, labels, loc='center left',
    fig.legend(handles, labels, loc="lower center",bbox_to_anchor=(0.5, -0.5),fancybox=True,ncol=len(grouped_columns), bbox_transform=axes[-2].transAxes) 

    # Save as one PDF
    fig.savefig(output_file, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved combined figure: {output_file}")

def plot_timeline(input_file, output_file):
    df = pd.read_csv(input_file,header=None)  
    cmap_name = df.iloc[0, 0]  # first column
    vmin = float(df.iloc[0, 1])  # second column
    vmax = float(df.iloc[0, 2])  # third column
    fig, ax = plt.subplots(figsize=(4, 0.4),constrained_layout=True)  # wide and short for horizontal colorbar
    # fig.subplots_adjust(right=0.5)  # adjust to make room for ticks

    # Create a ScalarMappable and hide axes
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = mpl.cm.ScalarMappable(cmap=cmap_name, norm=norm)
    sm.set_array([])

    # Add colorbar
    cbar = fig.colorbar(sm, cax=ax, orientation='horizontal',ticks=[vmin,0,vmax])
    # cbar.set_label(f"{cmap_name} ({vmin} to {vmax})")

    # Remove axes spines/ticks if you want it clean
    # ax.set_xticks([])
    # ax.set_yticks([])

    # Save as PDF
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    print(f"Saved colormap: {output_file}")


def add_outline(img_path, out_path="outlined.png", color=(0, 255, 0, 255), thickness=6):
    """
    Draws an outline (halo) outside the **outer shape** of an RGBA image,
    ignoring internal holes.
    
    Parameters:
        img_path (str): Path to the input image.
        out_path (str): Path to save the output image.
        color (tuple): RGBA color of the outline. Default green.
        thickness (int): Thickness of the outline (in pixels).
    """
    # Load image with alpha channel
    img = cv2.imread(img_path, cv2.IMREAD_UNCHANGED)
    if img is None:
        raise FileNotFoundError(f"Image not found: {img_path}")
    if img.shape[2] != 4:
        raise ValueError("Image must have an alpha channel (RGBA).")
    img = cv2.copyMakeBorder(img, 20, 20, 20, 20, cv2.BORDER_CONSTANT, value=[0,0,0,0])

    alpha = img[:, :, 3]
    _, mask = cv2.threshold(alpha, 0, 255, cv2.THRESH_BINARY)

    # Find only **external contours**
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Create empty outline image
    outline_img = np.zeros_like(img)

    # Draw contours with thickness
    cv2.drawContours(outline_img, contours, -1, color, thickness=thickness, lineType=cv2.LINE_AA)

    # Combine with original image: original pixels override outline
    alpha_mask = img[:, :, 3] > 0
    outline_img[alpha_mask] = img[alpha_mask]

    # Save result
    cv2.imwrite(out_path, outline_img)
    print(f"Outlined image saved to {out_path}")




def outline_file(args):
    figures_dir, f, color_map, name, case_dir = args
    directory, filename = os.path.split(f)  
    rel_dir = os.path.relpath(case_dir, start=os.path.dirname(case_dir))  
    out_dir = os.path.join(figures_dir, rel_dir)
    os.makedirs(out_dir, exist_ok=True)

    f_colored = os.path.join(out_dir, filename)
    rgba = tuple(int(255*x) for x in color_map[name])
    bgr_color = (rgba[2], rgba[1], rgba[0], rgba[3])  # (B,G,R,A)
    add_outline(f, f_colored, bgr_color, 16)

def postprocess_series(figures_dir, cases):
    """
    For each name and attribute, merge all case PNGs into one horizontal image
    and save as a PDF in colored_snapshots.
    """
    for name, _ in cases:
        # find all processed subfolders for this name
        case_dirs = sorted(glob.glob(os.path.join(figures_dir, f"{name}[0-9]*")))
        if not case_dirs:
            continue

        # Collect attributes by looking at first case folder
        first_case = case_dirs[0]
        attributes = [os.path.splitext(os.path.basename(f))[0]
                      for f in glob.glob(os.path.join(first_case, "*.png"))]

        for attr in attributes:
            # Collect images for this attribute across cases
            images = []
            for case_dir in case_dirs:
                f = os.path.join(case_dir, f"{attr}.png")
                if os.path.exists(f):
                    img = Image.open(f).convert("RGBA")
                    images.append(img)

            if not images:
                continue

            # Concatenate horizontally
            widths, heights = zip(*(img.size for img in images))
            total_width = sum(widths)
            max_height = max(heights)

            combined = Image.new("RGBA", (total_width, max_height))
            x_offset = 0
            for img in images:
                combined.paste(img, (x_offset, 0))
                x_offset += img.width

            background = Image.new("RGB", combined.size, (255, 255, 255))
            background.paste(combined, mask=combined.split()[3]) # use alpha channel as mask

            # Save as PDF
            out_pdf = os.path.join(figures_dir, f"{name}-{attr}.pdf") 
            background.save(out_pdf, "PDF", resolution=100.0) 
            print(f"Saved {out_pdf}")


def outline_pics(input_dir, cases):
    snapshots_dir = os.path.join(input_dir, "snapshots")
    figures_dir = os.path.join(input_dir, "colored_snapshots")
    os.makedirs(snapshots_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    color_map = assign_colors([name for name, _ in cases])

    for name, _ in cases:
        print(f"Outlining {name}")
        # find all subfolders like name[0-9]
        case_dirs = sorted(glob.glob(os.path.join(snapshots_dir, f"{name}[0-9]*")))
        for case_dir in case_dirs:
            files = sorted(glob.glob(os.path.join(case_dir, "*.png")))
            # for f in files:
            #     outline_file((figures_dir, f, color_map, name, case_dir))
            with ProcessPoolExecutor() as executor:
                executor.map(
                    outline_file,
                    [(figures_dir, f, color_map, name, case_dir) for f in files]
                )
        print("finished")
    postprocess_series(figures_dir, cases)


def get_config(path):
    name = os.path.basename(os.path.normpath(path))

    # ncols, figsize
    mapping = {
        "run_tc1": (3, (4, 2.1)),
        "run_tc2": (3, (4, 2.6)),
        "run_tc3.1": (4, (3, 2.5)),
        "run_tc3.2": (4, (3, 2.5)),
        "run_tc4": (3, (4, 2.3)),
        "run_tc5": (3, (4, 2)),
        "run_tc6": (3, (4, 2)),
        "run_tc7": (3, (4, 2)),
        "run_tc8": (2, (4, 4)),
        "run_tc9": (3, (4, 2)),
        "run_tc13": (3, (4, 4)),
        "run_tc14": (3, (4, 4)),
    }
    if name not in mapping:
        raise ValueError(f"Unsupported directory name: {name}")

    return mapping.get(name)



def plot_dir(input_dir, cases):
    data_dir = os.path.join(input_dir,"data")
    plot_dir = os.path.join(input_dir,"plots")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    dataframes = []
    for name, label in cases:
        print(f"Reading data of {name}")
        filepath = os.path.join(data_dir, name+"-measurements.csv")
        df = pd.read_csv(filepath)
        dataframes.append((name,label,df))


    if os.path.basename(os.path.normpath(input_dir)) == "run_tc5":
        plot_grouped_columns(dataframes,plot_dir)
    plot_bigpage(dataframes, plot_dir+"/all_plots.pdf",get_config(input_dir))
    plot_timeline(data_dir+"/colormap.csv",plot_dir+"/colormap.pdf")

def process_dir(input_dir):
    print(f"processing {input_dir}")
    df = pd.read_csv(input_dir + "/config.csv")
    df.columns = df.columns.str.strip()  # remove any whitespace in column names
    cases = [(row.name, row.label) for row in df.itertuples(index=False)]
    plot_dir(input_dir, cases)
    outline_pics(input_dir,cases)

if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = current_dir

    # find all run_tc* directories
    run_dirs = sorted(glob.glob(os.path.join(results_dir, "run_tc*")))

    # process_dir(run_dirs[0])
    if not run_dirs:
        print("No run_tc* directories found.")
    else:
        print(f"Found {len(run_dirs)} run directories")

        # run them in parallel
        with ProcessPoolExecutor() as executor:
            executor.map(process_dir, run_dirs)
