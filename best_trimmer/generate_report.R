# Function to safely load required packages
load_required_packages <- function() {
    required_packages <- c(
        "ggplot2",     # Core plotting
        "dplyr",       # Data manipulation
        "tidyr",       # Data reshaping
        "ggridges",    # Ridgeline plots
        "hexbin",      # Hexagonal binning
        "RColorBrewer", # Color palettes
        "ggtext",      # Rich text support
        "MetBrewer",   # Metropolitan color palettes
        "ggdist",      # Distribution visualization
        "patchwork"    # For inset_element functionality
    )
    
    # Check and install missing packages
    missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
    if(length(missing_packages) > 0) {
        warning("Installing missing packages: ", paste(missing_packages, collapse=", "))
        install.packages(missing_packages)
    }
    
    # Safely load each package
    for(pkg in required_packages) {
        tryCatch(
            suppressPackageStartupMessages(library(pkg, character.only = TRUE)),
            error = function(e) {
                stop(paste("Failed to load package:", pkg, "\nError:", e$message))
            }
        )
    }
}

# Load packages
load_required_packages()

# Memory management function
manage_memory <- function() {
    gc()
}

# Read the data with safeguard
tryCatch({
    data <- read.csv("read_metrics.csv")
}, error=function(e) {
    stop("Error reading data: ", e$message)
})

# Process in chunks
create_and_save_plots <- function(data) {
    # Style settings
    bg_color <- "grey97"
    font_family <- "Arial"  # Changed from "Fira Sans" to "Arial"
    
    # Create list to store plots
    plot_list <- list()
    
    # Sequence length plot with both density and statistical overlays
    plot_list[[1]] <- ggplot(data, aes(x=Sequence.Length, y=Tool)) +
        geom_density_ridges(aes(fill = Tool), alpha = 0.3, scale = 0.7) +
        stat_halfeye(aes(fill = Tool), fill_type = "segments", alpha = 0.3, scale = 0.7) +
        stat_interval(aes(color = Tool)) +
        stat_summary(geom = "point", fun = median) +
        geom_hline(yintercept = median(data$Sequence.Length), col = "grey30", lty = "dashed") +
        scale_y_discrete(labels = toupper) +
        scale_color_manual(values = met.brewer("VanGogh3")) +
        scale_fill_manual(values = met.brewer("VanGogh3")) +
        theme_minimal(base_family = font_family) +
        theme(
            plot.background = element_rect(color = NA, fill = bg_color),
            panel.grid = element_blank(),
            panel.grid.major.x = element_line(linewidth = 0.1, color = "grey75"),
            plot.title = element_text(family = font_family, face = "bold"),  # Changed from Fira Sans SemiBold
            plot.title.position = "plot",
            axis.text.y = element_text(hjust = 0, family = font_family, face = "bold"),  # Changed from Fira Sans SemiBold
            legend.position = "none"
        ) +
        labs(title = "Sequence Length Distribution by Tool",
             x = "Sequence Length", y = NULL)

    # GC content ridgeline plot
    plot_list[[2]] <- ggplot(data, aes(x=GC.Content, y=Tool, fill=Tool)) +
        geom_density_ridges(scale=0.9, alpha=0.6, bandwidth=0.00587) +
        theme_ridges() +
        theme(legend.position="none") +
        labs(title="GC Content Distribution by Tool",
             x="GC Content", y="Tool")

    # Base composition ridgeline plot
    base_data <- data %>%
        select(Tool, starts_with("Base_")) %>%
        gather(key="Base", value="Percentage", -Tool)
    
    plot_list[[3]] <- ggplot(base_data, aes(x=Percentage, y=Tool, fill=Base)) +
        geom_density_ridges(scale=0.9, alpha=0.6, bandwidth=0.346) +
        theme_ridges() +
        labs(title="Base Composition Distribution by Tool",
             x="Percentage", y="Tool")

    # Hexbin plot for entropy vs read length
    plot_list[[4]] <- ggplot(data, aes(x=Sequence.Length, y=Entropy)) +
        geom_hex(bins=40) +
        facet_wrap(~Tool) +
        scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
        theme_minimal() +
        labs(title="Sequence Complexity (Entropy) vs Read Length by Trimmer",
             x="Sequence Length",
             y="Shannon Entropy") +
        theme(legend.position="right")

    # Create and add legend plot
    legend_data <- subset(data, Tool == levels(data$Tool)[1])
    p_legend <- ggplot(legend_data, aes(x=Tool, y=Sequence.Length)) +
        stat_halfeye(fill_type = "segments", alpha = 0.3) +
        stat_interval() +
        stat_summary(geom = "point", fun = median) +
        annotate(
            "richtext",
            x = c(0.8, 0.8, 0.8, 1.4, 1.8),
            y = c(50, 150, 100, 75, 125),
            label = c("50% of lengths\nfall within this range", 
                     "95% of lengths", "80% of lengths", 
                     "Median", "Distribution\nof lengths"),
            fill = NA, label.size = 0, family = font_family, size = 2, vjust = 1
        ) +
        theme_void(base_family = font_family) +
        theme(plot.background = element_rect(color = "grey30", 
                                           linewidth = 0.2,  # Changed from size to linewidth
                                           fill = bg_color))

    # Save plots with legend
    pdf("sequence_report.pdf", width=12, height=20)
    for(i in seq_along(plot_list)) {
        if(i == 1) {
            print(plot_list[[i]] + 
                  inset_element(p_legend, l = 0.6, r = 1.0, t = 0.99, b = 0.7))
        } else {
            print(plot_list[[i]])
        }
    }
    dev.off()
    
    # Clean up
    rm(plot_list, base_data, legend_data)
    manage_memory()
}

# Main execution
create_and_save_plots(data)
manage_memory()
