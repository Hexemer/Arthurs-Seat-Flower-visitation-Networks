##################################
#Plant-pollinator Network Figures
##################################

# Note: objects labelled "pollinator" refer to observed flower visitors.
# Pollination effectiveness was not directly measured.

#Import data

pollinators <- read.csv(
  "Raw data/Pollinators.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Date",
  "Species",
  "Genus",
  "Family",
  "Order",
  "Visits",
  "Block",
  "Replication",
  "Time.Period",
  "Site",
  "Plant.ID",
  "Number.of.Flowers",
  "Abundance",
  "Start.Time"
)


missing_columns <- setdiff(
  required_columns,
  names(pollinators)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "These required columns are missing:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


#Plant names and categories

plant_key <- c(
  "1"  = "Yellow flatweed",
  "2"  = "Maiden pink",
  "3"  = "Thyme",
  "4"  = "Bird's-foot trefoil",
  "5"  = "Sticky catchfly",
  "6"  = "Heath bedstraw",
  "7"  = "Yellow flatweed",
  "8"  = "Thyme",
  "9"  = "Bird's-foot trefoil",
  "10" = "Sticky catchfly",
  "11" = "Woodland germander",
  "12" = "Lady's bedstraw",
  "13" = "Lady's bedstraw",
  "14" = "Plume thistle"
)


# Restored plants:
# Plant 2 = Maiden pink
# Plants 5 and 10 = Sticky catchfly

restored_ids <- c(
  "2",
  "5",
  "10"
)


restored_colour <- "#E7298A"
other_plant_colour <- "black"

plant_order <- as.character(
  1:14
)


#Clean Text

clean_text <- function(x) {
  
  x <- as.character(x)
  
  x <- gsub(
    "\u00A0",
    " ",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "[[:space:]]+",
    " ",
    x
  )
  
  trimws(x)
}


#Prepare Network Data

network_data <- pollinators %>%
  
  mutate(
    
    Species = clean_text(Species),
    Genus = clean_text(Genus),
    Family = clean_text(Family),
    Order = clean_text(Order),
    
    # Standardise family names
    Family = case_when(
      tolower(Family) == "syrphidae" ~ "Syrphidae",
      TRUE ~ Family
    ),
    
    # Use the most precise available pollinator identity
    Pollinator = case_when(
      
      !is.na(Species) &
        Species != "" ~
        Species,
      
      !is.na(Genus) &
        Genus != "" ~
        paste0(
          Genus,
          " sp."
        ),
      
      !is.na(Family) &
        Family != "" ~
        paste0(
          Family,
          " family"
        ),
      
      !is.na(Order) &
        Order != "" ~
        paste0(
          Order,
          " order"
        ),
      
      TRUE ~ NA_character_
    ),
    
    # Family used for colours
    Family.Network = case_when(
      
      !is.na(Family) &
        Family != "" ~
        Family,
      
      !is.na(Order) &
        Order != "" ~
        paste0(
          Order,
          " unidentified"
        ),
      
      TRUE ~ "Unidentified"
    ),
    
    Pollinator = clean_text(Pollinator),
    Family.Network = clean_text(Family.Network),
    
    Family.Network = case_when(
      tolower(Family.Network) == "syrphidae" ~
        "Syrphidae",
      TRUE ~ Family.Network
    ),
    
    Plant.ID = as.character(Plant.ID),
    Site = as.character(Site),
    Block = as.character(Block),
    Replication = as.character(Replication),
    Time.Period = as.character(Time.Period),
    
    Plant.Name = unname(
      plant_key[Plant.ID]
    ),
    
    Visit.Rate = case_when(
      
      !is.na(Number.of.Flowers) &
        Number.of.Flowers > 0 &
        !is.na(Visits) ~
        Visits / Number.of.Flowers,
      
      TRUE ~ NA_real_
    )
  )


#Interaction and floral data

interaction_data <- network_data %>%
  
  filter(
    Visits > 0,
    !is.na(Visit.Rate),
    Visit.Rate > 0,
    !is.na(Pollinator),
    Pollinator != "",
    !is.na(Plant.ID),
    Plant.ID != "",
    !is.na(Site),
    Site != "",
    !is.na(Block),
    Block != "",
    !is.na(Time.Period),
    Time.Period != ""
  )


floral_data <- network_data %>%
  
  filter(
    !is.na(Plant.ID),
    Plant.ID != "",
    !is.na(Abundance),
    Abundance > 0,
    !is.na(Site),
    Site != "",
    !is.na(Block),
    Block != "",
    !is.na(Time.Period),
    Time.Period != ""
  ) %>%
  
  distinct(
    Date,
    Site,
    Block,
    Replication,
    Time.Period,
    Plant.ID,
    Start.Time,
    .keep_all = TRUE
  )


#Time-period order

time_period_order <- network_data %>%
  
  filter(
    !is.na(Time.Period),
    Time.Period != ""
  ) %>%
  
  distinct(Time.Period) %>%
  
  pull(Time.Period) %>%
  
  as.character()


numeric_time_periods <- suppressWarnings(
  as.numeric(time_period_order)
)


if (all(!is.na(numeric_time_periods))) {
  
  time_period_order <- time_period_order[
    order(numeric_time_periods)
  ]
}


if (length(time_period_order) != 3) {
  
  warning(
    paste(
      "Expected three time periods but found",
      length(time_period_order),
      ":",
      paste(
        time_period_order,
        collapse = ", "
      )
    )
  )
}


#Family colours

families <- sort(
  unique(
    interaction_data$Family.Network
  )
)


family_colours <- grDevices::hcl.colors(
  n = max(
    length(families),
    1
  ),
  palette = "Dark 3"
)


family_colours <- family_colours[
  seq_along(families)
]


names(family_colours) <- families


#Abbreviate pollinator names

abbreviate_pollinator <- function(x) {
  
  x <- as.character(x)
  
  vapply(
    
    x,
    
    function(label) {
      
      # Keep unidentified and higher-level labels unchanged
      if (
        grepl(
          " family$| order$| unidentified$| sp\\.$",
          label,
          ignore.case = TRUE
        )
      ) {
        
        return(label)
      }
      
      words <- strsplit(
        label,
        "[[:space:]]+"
      )[[1]]
      
      # Abbreviate normal binomial names
      if (length(words) >= 2) {
        
        return(
          paste0(
            substr(
              words[1],
              1,
              1
            ),
            ". ",
            paste(
              words[-1],
              collapse = " "
            )
          )
        )
      }
      
      label
    },
    
    character(1)
  )
}


#Floral visitors

#Floral visitors are ordered from greatest to smallest
# total interaction strength within each block.

get_pollinator_order <- function(
    block_interactions
) {
  
  block_interactions %>%
    
    group_by(
      Pollinator
    ) %>%
    
    summarise(
      
      Total.Rate = sum(
        Visit.Rate,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    
    arrange(
      desc(Total.Rate),
      Pollinator
    ) %>%
    
    pull(Pollinator)
}


# ==========================================================
# 10. BUILD FIXED BLOCK LAYOUT
# ==========================================================

build_block_layout <- function(
    block_interactions,
    block_flowers,
    pollinator_order
) {
  
  
  # --------------------------------------------------------
  # Plant widths across the block
  # --------------------------------------------------------
  
  block_plant_totals <- block_flowers %>%
    
    group_by(
      Plant.ID
    ) %>%
    
    summarise(
      
      Floral.Abundance = mean(
        Abundance,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    
    filter(
      is.finite(Floral.Abundance),
      Floral.Abundance > 0
    )
  
  
  linked_plants <- unique(
    block_interactions$Plant.ID
  )
  
  
  present_plants <- union(
    block_plant_totals$Plant.ID,
    linked_plants
  )
  
  
  present_plants <- plant_order[
    plant_order %in% present_plants
  ]
  
  
  block_plant_totals <- tibble(
    Plant.ID = present_plants
  ) %>%
    
    left_join(
      block_plant_totals,
      by = "Plant.ID"
    ) %>%
    
    mutate(
      
      Floral.Abundance = replace_na(
        Floral.Abundance,
        0
      )
    )
  
  
  # --------------------------------------------------------
  # Pollinator widths across the block
  # --------------------------------------------------------
  
  block_pollinator_totals <- block_interactions %>%
    
    group_by(
      Pollinator
    ) %>%
    
    summarise(
      
      Pollinator.Total = sum(
        Visit.Rate,
        na.rm = TRUE
      ),
      
      Family.Network = first(
        Family.Network
      ),
      
      .groups = "drop"
    )
  
  
  present_pollinators <- pollinator_order[
    pollinator_order %in%
      block_pollinator_totals$Pollinator
  ]
  
  
  block_pollinator_totals <- tibble(
    Pollinator = present_pollinators
  ) %>%
    
    left_join(
      block_pollinator_totals,
      by = "Pollinator"
    ) %>%
    
    filter(
      !is.na(Pollinator.Total),
      Pollinator.Total > 0
    )
  
  
  # --------------------------------------------------------
  # Transform bar widths
  # --------------------------------------------------------
  
  plant_widths <- sqrt(
    block_plant_totals$Floral.Abundance
  )
  
  
  pollinator_widths <- sqrt(
    block_pollinator_totals$Pollinator.Total
  )
  
  
  if (
    length(plant_widths) > 0 &&
    max(
      plant_widths,
      na.rm = TRUE
    ) > 0
  ) {
    
    plant_widths <- pmax(
      
      plant_widths,
      
      max(
        plant_widths,
        na.rm = TRUE
      ) * 0.025
    )
  }
  
  
  if (
    length(pollinator_widths) > 0 &&
    max(
      pollinator_widths,
      na.rm = TRUE
    ) > 0
  ) {
    
    pollinator_widths <- pmax(
      
      pollinator_widths,
      
      max(
        pollinator_widths,
        na.rm = TRUE
      ) * 0.07
    )
  }
  
  
  # --------------------------------------------------------
  # Original plant widths
  # --------------------------------------------------------
  
  plant_gap <- max(
    sum(plant_widths) * 0.01,
    0.02
  )
  
  
  plant_starts <- cumsum(
    c(
      0,
      head(
        plant_widths,
        -1
      ) + plant_gap
    )
  )
  
  
  plant_ends <- plant_starts +
    plant_widths
  
  
  plant_centres <- (
    plant_starts +
      plant_ends
  ) / 2
  
  
  # --------------------------------------------------------
  # Fixed pollinator positions
  # --------------------------------------------------------
  
  pollinator_gap <- max(
    sum(pollinator_widths) * 0.008,
    0.015
  )
  
  
  pollinator_starts <- cumsum(
    c(
      0,
      head(
        pollinator_widths,
        -1
      ) + pollinator_gap
    )
  )
  
  
  pollinator_ends <- pollinator_starts +
    pollinator_widths
  
  
  pollinator_centres <- (
    pollinator_starts +
      pollinator_ends
  ) / 2
  
  
  # --------------------------------------------------------
  # Scale both levels to the panel width
  # --------------------------------------------------------
  
  plotting_width <- 1
  
  
  plant_total_width <- if (
    length(plant_ends) > 0
  ) {
    
    max(
      plant_ends,
      na.rm = TRUE
    )
    
  } else {
    
    1
  }
  
  
  pollinator_total_width <- if (
    length(pollinator_ends) > 0
  ) {
    
    max(
      pollinator_ends,
      na.rm = TRUE
    )
    
  } else {
    
    1
  }
  
  
  plant_scale <-
    plotting_width /
    plant_total_width
  
  
  pollinator_scale <-
    plotting_width /
    pollinator_total_width
  
  
  plant_starts <-
    plant_starts *
    plant_scale
  
  
  plant_ends <-
    plant_ends *
    plant_scale
  
  
  plant_centres <-
    plant_centres *
    plant_scale
  
  
  pollinator_starts <-
    pollinator_starts *
    pollinator_scale
  
  
  pollinator_ends <-
    pollinator_ends *
    pollinator_scale
  
  
  pollinator_centres <-
    pollinator_centres *
    pollinator_scale
  
  
  plant_positions <- block_plant_totals %>%
    
    mutate(
      Start = plant_starts,
      End = plant_ends,
      Centre = plant_centres
    )
  
  
  pollinator_positions <- block_pollinator_totals %>%
    
    mutate(
      Start = pollinator_starts,
      End = pollinator_ends,
      Centre = pollinator_centres
    )
  
  
  # --------------------------------------------------------
  # Consistent link-width scale within each block
  # --------------------------------------------------------
  
  overall_links <- block_interactions %>%
    
    group_by(
      Site,
      Plant.ID,
      Pollinator
    ) %>%
    
    summarise(
      
      Weight = sum(
        Visit.Rate,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  
  period_links <- block_interactions %>%
    
    group_by(
      Site,
      Time.Period,
      Plant.ID,
      Pollinator
    ) %>%
    
    summarise(
      
      Weight = sum(
        Visit.Rate,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  
  all_panel_weights <- c(
    overall_links$Weight,
    period_links$Weight
  )
  
  
  all_panel_weights <- all_panel_weights[
    is.finite(all_panel_weights) &
      all_panel_weights > 0
  ]
  
  
  transformed_weights <- sqrt(
    all_panel_weights
  )
  
  
  if (length(transformed_weights) == 0) {
    
    link_min <- 0
    link_max <- 0
    
  } else {
    
    link_min <- min(
      transformed_weights
    )
    
    link_max <- max(
      transformed_weights
    )
  }
  
  
  list(
    
    plotting_width = plotting_width,
    
    plant_positions =
      plant_positions,
    
    pollinator_positions =
      pollinator_positions,
    
    link_min =
      link_min,
    
    link_max =
      link_max
  )
}


# ==========================================================
# 11. PACK PLANTS PRESENT IN ONE PANEL
# ==========================================================

pack_panel_plants <- function(
    flowers,
    interactions,
    block_layout
) {
  
  panel_plants <- union(
    unique(flowers$Plant.ID),
    unique(interactions$Plant.ID)
  )
  
  
  panel_plants <- plant_order[
    plant_order %in% panel_plants
  ]
  
  
  panel_positions <-
    block_layout$plant_positions %>%
    
    filter(
      Plant.ID %in% panel_plants
    ) %>%
    
    mutate(
      Panel.Width = End - Start
    )
  
  
  if (nrow(panel_positions) == 0) {
    
    return(panel_positions)
  }
  
  
  raw_widths <- panel_positions$Panel.Width
  
  
  panel_gap <- max(
    sum(raw_widths) * 0.018,
    0.012
  )
  
  
  new_starts <- cumsum(
    c(
      0,
      head(
        raw_widths,
        -1
      ) + panel_gap
    )
  )
  
  
  new_ends <- new_starts +
    raw_widths
  
  
  total_width <- max(
    new_ends,
    na.rm = TRUE
  )
  
  
  if (
    !is.finite(total_width) ||
    total_width <= 0
  ) {
    
    total_width <- 1
  }
  
  
  panel_scale <-
    block_layout$plotting_width /
    total_width
  
  
  panel_positions %>%
    
    mutate(
      
      Start =
        new_starts *
        panel_scale,
      
      End =
        new_ends *
        panel_scale,
      
      Centre = (
        Start + End
      ) / 2
    )
}


# ==========================================================
# 12. DRAW ONE NETWORK PANEL
# ==========================================================

plot_custom_network <- function(
    interactions,
    flowers,
    panel_title,
    block_layout,
    show_pollinator_labels = FALSE
) {
  
  
  plot(
    NA,
    xlim = c(
      0,
      block_layout$plotting_width
    ),
    ylim = c(
      0,
      1.35
    ),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i"
  )
  
  
  title(
    main = panel_title,
    adj = 0,
    line = 1,
    cex.main = 1.35,
    font.main = 2
  )
  
  
  # --------------------------------------------------------
  # Empty panel
  # --------------------------------------------------------
  
  if (
    nrow(interactions) == 0 ||
    sum(
      interactions$Visit.Rate,
      na.rm = TRUE
    ) <= 0
  ) {
    
    text(
      x = block_layout$plotting_width / 2,
      y = 0.5,
      labels = "No flower-visitor visits recorded",
      cex = 0.9
    )
    
    return(
      invisible(NULL)
    )
  }
  
  
  # --------------------------------------------------------
  # Panel-specific plant positions
  # --------------------------------------------------------
  
  panel_plant_positions <- pack_panel_plants(
    flowers = flowers,
    interactions = interactions,
    block_layout = block_layout
  )
  
  
  # --------------------------------------------------------
  # Create interaction links
  # --------------------------------------------------------
  
  links <- interactions %>%
    
    group_by(
      Plant.ID,
      Pollinator,
      Family.Network
    ) %>%
    
    summarise(
      
      Weight = sum(
        Visit.Rate,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    
    left_join(
      
      panel_plant_positions %>%
        
        select(
          Plant.ID,
          Plant.X = Centre
        ),
      
      by = "Plant.ID"
    ) %>%
    
    left_join(
      
      block_layout$pollinator_positions %>%
        
        select(
          Pollinator,
          Pollinator.X = Centre
        ),
      
      by = "Pollinator"
    ) %>%
    
    filter(
      !is.na(Plant.X),
      !is.na(Pollinator.X)
    )
  
  
  # --------------------------------------------------------
  # Calculate link widths
  # --------------------------------------------------------
  
  transformed_weights <- sqrt(
    links$Weight
  )
  
  
  if (
    block_layout$link_max <=
    block_layout$link_min
  ) {
    
    line_widths <- rep(
      2.5,
      length(transformed_weights)
    )
    
  } else {
    
    line_widths <- 1.1 +
      3.4 * (
        
        (
          transformed_weights -
            block_layout$link_min
        ) /
          
          (
            block_layout$link_max -
              block_layout$link_min
          )
      )
    
    
    line_widths <- pmax(
      1.1,
      pmin(
        4.5,
        line_widths
      )
    )
  }
  
  
  # --------------------------------------------------------
  # Draw links
  # --------------------------------------------------------
  
  for (i in seq_len(nrow(links))) {
    
    link_colour <- family_colours[
      links$Family.Network[i]
    ]
    
    
    if (
      length(link_colour) == 0 ||
      is.na(link_colour)
    ) {
      
      link_colour <- "grey50"
    }
    
    
    segments(
      
      x0 = links$Plant.X[i],
      y0 = 0.15,
      
      x1 = links$Pollinator.X[i],
      y1 = 0.82,
      
      col = grDevices::adjustcolor(
        link_colour,
        alpha.f = 0.40
      ),
      
      lwd = line_widths[i],
      lend = "butt"
    )
  }
  
  
  # --------------------------------------------------------
  # Plant bars
  # --------------------------------------------------------
  
  plant_bars <- panel_plant_positions
  
  
  plant_bar_colours <- ifelse(
    plant_bars$Plant.ID %in%
      restored_ids,
    restored_colour,
    other_plant_colour
  )
  
  
  rect(
    xleft = plant_bars$Start,
    ybottom = 0.08,
    xright = plant_bars$End,
    ytop = 0.15,
    col = plant_bar_colours,
    border = plant_bar_colours
  )
  
  
  # --------------------------------------------------------
  # Plant names
  # --------------------------------------------------------
  
  plant_labels <- unname(
    plant_key[
      plant_bars$Plant.ID
    ]
  )
  
  
  text(
    x = plant_bars$Centre,
    y = 0.065,
    labels = plant_labels,
    srt = 45,
    adj = c(
      1,
      0.5
    ),
    cex = 0.68,
    font = 2,
    xpd = NA
  )
  
  
  # --------------------------------------------------------
  # Pollinator bars
  # --------------------------------------------------------
  
  panel_pollinators <- unique(
    links$Pollinator
  )
  
  
  pollinator_bars <-
    block_layout$pollinator_positions %>%
    
    filter(
      Pollinator %in%
        panel_pollinators
    )
  
  
  pollinator_bar_colours <-
    family_colours[
      pollinator_bars$Family.Network
    ]
  
  
  pollinator_bar_colours[
    is.na(pollinator_bar_colours)
  ] <- "grey50"
  
  
  rect(
    xleft = pollinator_bars$Start,
    ybottom = 0.82,
    xright = pollinator_bars$End,
    ytop = 0.93,
    col = pollinator_bar_colours,
    border = pollinator_bar_colours
  )
  
  
  # --------------------------------------------------------
  # Pollinator names on overall panels only
  # --------------------------------------------------------
  
  if (show_pollinator_labels) {
    
    pollinator_labels <- abbreviate_pollinator(
      pollinator_bars$Pollinator
    )
    
    
    text(
      x = pollinator_bars$Centre,
      y = 0.945,
      labels = pollinator_labels,
      srt = 90,
      adj = c(
        0,
        0.5
      ),
      cex = 0.72,
      col = pollinator_bar_colours,
      xpd = NA
    )
  }
  
  
  invisible(NULL)
}


# ==========================================================
# 13. LEGEND
# ==========================================================

plot_legend_panel <- function(
    block_interactions
) {
  
  
  plot.new()
  
  
  plot.window(
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    )
  )
  
  
  used_families <- sort(
    unique(
      block_interactions$Family.Network
    )
  )
  
  
  used_colours <- family_colours[
    used_families
  ]
  
  
  text(
    x = 0.5,
    y = 0.90,
    labels = "Flower-visitor family",
    font = 2,
    cex = 1.4
  )
  
  
  if (length(used_families) > 0) {
    
    number_of_columns <- min(
      4,
      length(used_families)
    )
    
    
    number_of_rows <- ceiling(
      length(used_families) /
        number_of_columns
    )
    
    
    legend_x <- rep(
      
      seq(
        0.08,
        0.82,
        length.out =
          number_of_columns
      ),
      
      times =
        number_of_rows
      
    )[seq_along(used_families)]
    
    
    legend_y <- rep(
      
      seq(
        0.72,
        0.48,
        length.out =
          number_of_rows
      ),
      
      each =
        number_of_columns
      
    )[seq_along(used_families)]
    
    
    rect(
      xleft = legend_x,
      ybottom = legend_y - 0.035,
      xright = legend_x + 0.025,
      ytop = legend_y + 0.035,
      col = used_colours,
      border = used_colours
    )
    
    
    text(
      x = legend_x + 0.035,
      y = legend_y,
      labels = used_families,
      adj = c(
        0,
        0.5
      ),
      cex = 0.95
    )
  }
  
  
  text(
    x = 0.5,
    y = 0.30,
    labels = "Plant category",
    font = 2,
    cex = 1.4
  )
  
  
  plant_x <- c(
    0.27,
    0.61
  )
  
  
  rect(
    xleft = plant_x,
    ybottom = 0.12,
    xright = plant_x + 0.025,
    ytop = 0.19,
    col = c(
      other_plant_colour,
      restored_colour
    ),
    border = c(
      other_plant_colour,
      restored_colour
    )
  )
  
  
  text(
    x = plant_x + 0.035,
    y = 0.155,
    labels = c(
      "Non-restored plants",
      "Restored plants"
    ),
    adj = c(
      0,
      0.5
    ),
    cex = 1
  )
  
  
  invisible(NULL)
}


# ==========================================================
# 14. CREATE COMPLETE BLOCK FIGURE
# ==========================================================

plot_complete_block <- function(
    block_to_plot
) {
  
  
  block_to_plot <- as.character(
    block_to_plot
  )
  
  
  block_interactions <-
    interaction_data %>%
    
    filter(
      Block == block_to_plot
    )
  
  
  block_flowers <-
    floral_data %>%
    
    filter(
      Block == block_to_plot
    )
  
  
  if (nrow(block_interactions) == 0) {
    
    stop(
      paste(
        "No interaction data found for Block",
        block_to_plot
      )
    )
  }
  
  
  block_pollinator_order <-
    get_pollinator_order(
      block_interactions
    )
  
  
  block_layout <-
    build_block_layout(
      
      block_interactions =
        block_interactions,
      
      block_flowers =
        block_flowers,
      
      pollinator_order =
        block_pollinator_order
    )
  
  
  layout_matrix <- matrix(
    
    c(
      1, 2,
      3, 4,
      5, 6,
      7, 8,
      9, 9
    ),
    
    byrow = TRUE,
    ncol = 2
  )
  
  
  layout(
    
    layout_matrix,
    
    heights = c(
      1.25,
      1,
      1,
      1,
      0.72
    )
  )
  
  
  # Larger upper and lower margins prevent labels being cut off
  
  par(
    
    mar = c(
      4.2,
      0.8,
      7.2,
      0.8
    ),
    
    oma = c(
      0.5,
      0.5,
      4.5,
      0.5
    ),
    
    xpd = NA
  )
  
  
  sites <- c(
    "1",
    "2"
  )
  
  
  # --------------------------------------------------------
  # Overall panels
  # --------------------------------------------------------
  
  for (site_to_plot in sites) {
    
    
    site_interactions <-
      block_interactions %>%
      
      filter(
        Site == site_to_plot
      )
    
    
    site_flowers <-
      block_flowers %>%
      
      filter(
        Site == site_to_plot
      )
    
    
    plot_custom_network(
      
      interactions =
        site_interactions,
      
      flowers =
        site_flowers,
      
      panel_title = paste(
        "Site",
        site_to_plot,
        "– Overall"
      ),
      
      block_layout =
        block_layout,
      
      show_pollinator_labels =
        TRUE
    )
  }
  
  
  # --------------------------------------------------------
  # Time-period panels
  # --------------------------------------------------------
  
  for (time_to_plot in time_period_order) {
    
    for (site_to_plot in sites) {
      
      
      panel_interactions <-
        block_interactions %>%
        
        filter(
          Site == site_to_plot,
          Time.Period == time_to_plot
        )
      
      
      panel_flowers <-
        block_flowers %>%
        
        filter(
          Site == site_to_plot,
          Time.Period == time_to_plot
        )
      
      
      plot_custom_network(
        
        interactions =
          panel_interactions,
        
        flowers =
          panel_flowers,
        
        panel_title = paste(
          "Site",
          site_to_plot,
          "– Time period",
          time_to_plot
        ),
        
        block_layout =
          block_layout,
        
        show_pollinator_labels =
          FALSE
      )
    }
  }
  
  
  # --------------------------------------------------------
  # Legend panel
  # --------------------------------------------------------
  
  par(
    mar = c(
      0.3,
      0.3,
      0.3,
      0.3
    )
  )
  
  
  plot_legend_panel(
    block_interactions
  )
  
  
  # --------------------------------------------------------
  # Overall block title
  # --------------------------------------------------------
  
  mtext(
    paste(
      "Block",
      block_to_plot
    ),
    side = 3,
    outer = TRUE,
    line = 0.6,
    cex = 2.3,
    font = 2
  )
}


# ==========================================================
# 15. CLOSE OLD GRAPHICS DEVICES
# ==========================================================

while (dev.cur() > 1) {
  
  dev.off()
}


# ==========================================================
# 16. PREVIEW BLOCK 1
# ==========================================================

dev.new(
  width = 12,
  height = 15
)

plot_complete_block(
  "1"
)


# ==========================================================
# 17. PREVIEW BLOCK 2
# ==========================================================

dev.new(
  width = 12,
  height = 15
)

plot_complete_block(
  "2"
)


# ==========================================================
# 18. SAVE BLOCK 1
# ==========================================================

pdf(
  "Block_1_network_updated.pdf",
  width = 12,
  height = 15,
  family = "sans"
)

plot_complete_block(
  "1"
)

dev.off()


# ==========================================================
# 19. SAVE BLOCK 2
# ==========================================================

pdf(
  "Block_2_network_updated.pdf",
  width = 12,
  height = 15,
  family = "sans"
)

plot_complete_block(
  "2"
)

dev.off()

