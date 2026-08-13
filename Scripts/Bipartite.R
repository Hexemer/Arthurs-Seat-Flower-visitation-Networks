##################################
#Plant-flower visitor Network Figures
##################################

library(dplyr)
library(tidyr)
library(tibble)


#Plant species and categories

plant_key <- c(
  "1"  = "Hypochaeris radicata",
  "2"  = "Dianthus deltoides",
  "3"  = "Thymus drucei",
  "4"  = "Lotus corniculatus",
  "5"  = "Silene viscaria",
  "6"  = "Galium saxatile",
  "7"  = "Hypochaeris radicata",
  "8"  = "Thymus drucei",
  "9"  = "Lotus corniculatus",
  "10" = "Silene viscaria",
  "11" = "Teucrium scorodonia",
  "12" = "Galium verum",
  "13" = "Galium verum",
  "14" = "Cirsium arvense"
)

plant_order <- c(
  "Hypochaeris radicata",
  "Dianthus deltoides",
  "Thymus drucei",
  "Lotus corniculatus",
  "Silene viscaria",
  "Galium saxatile",
  "Teucrium scorodonia",
  "Galium verum",
  "Cirsium arvense"
)

restored_plants <- c(
  "Dianthus deltoides",
  "Silene viscaria"
)

restored_colour <- "#E7298A"

other_plant_colour <- "black"

time_period_order <- c(
  "1",
  "2",
  "3"
)


#Clean network text

clean_network_text <- function(x) {
  
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


#Prepare network data

network_figure_data <- pollinators %>%
  mutate(
    
    Species = clean_network_text(Species),
    Genus = clean_network_text(Genus),
    Family = clean_network_text(Family),
    Order = clean_network_text(Order),
    Group = clean_network_text(Group),
    
    Plant.ID = as.character(
      Plant.ID
    ),
    
    Site = as.character(
      Site
    ),
    
    Time.Period = as.character(
      Time.Period
    ),
    
    Plant.Species = unname(
      plant_key[
        Plant.ID
      ]
    ),
    
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
      
      TRUE ~
        "Unidentified"
    ),
    
    #Internal name remains Pollinator but represents flower visitors
    
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
      
      TRUE ~
        NA_character_
    ),
    
    Visit.Rate = case_when(
      
      !is.na(Visits) &
        !is.na(Number.of.Flowers) &
        Number.of.Flowers > 0 ~
        Visits / Number.of.Flowers,
      
      TRUE ~
        NA_real_
    )
  )


#Interaction data

network_interactions <- network_figure_data %>%
  filter(
    Group != "None",
    Visits > 0,
    !is.na(Visit.Rate),
    Visit.Rate > 0,
    !is.na(Plant.Species),
    Plant.Species != "",
    !is.na(Pollinator),
    Pollinator != "",
    !is.na(Site),
    Site != "",
    !is.na(Time.Period),
    Time.Period != ""
  )


#Floral data

network_flowers <- network_figure_data %>%
  filter(
    !is.na(Plant.Species),
    Plant.Species != "",
    !is.na(Abundance),
    Abundance > 0,
    !is.na(Site),
    Site != "",
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


#Flower visitor colours

pollinator_families <- sort(
  unique(
    network_interactions$Family.Network
  )
)

family_colours <- grDevices::hcl.colors(
  n = length(
    pollinator_families
  ),
  palette = "Dark 3"
)

names(
  family_colours
) <- pollinator_families


#Labels

abbreviate_pollinator <- function(x) {
  
  vapply(
    as.character(x),
    
    function(label) {
      
      if (
        grepl(
          " family$| order$| unidentified$| sp\\.$",
          label,
          ignore.case = TRUE
        )
      ) {
        
        return(
          label
        )
      }
      
      words <- strsplit(
        label,
        "[[:space:]]+"
      )[[1]]
      
      if (
        length(words) >= 2
      ) {
        
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


abbreviate_plant <- function(x) {
  
  vapply(
    as.character(x),
    
    function(label) {
      
      words <- strsplit(
        label,
        "[[:space:]]+"
      )[[1]]
      
      if (
        length(words) >= 2
      ) {
        
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


#Bar positions

make_bar_positions <- function(
    labels,
    weights,
    left = 0.04,
    right = 0.88,
    gap = 0.012
) {
  
  if (
    length(labels) == 0
  ) {
    
    return(
      tibble(
        Label = character(),
        Start = numeric(),
        End = numeric(),
        Centre = numeric()
      )
    )
  }
  
  weights[
    !is.finite(weights) |
      weights <= 0
  ] <- 0
  
  transformed <- sqrt(
    weights
  )
  
  if (
    max(
      transformed,
      na.rm = TRUE
    ) == 0
  ) {
    
    transformed <- rep(
      1,
      length(labels)
    )
  }
  
  minimum_width <- max(
    transformed
  ) * 0.05
  
  transformed <- pmax(
    transformed,
    minimum_width
  )
  
  available_width <-
    right -
    left -
    gap * (
      length(labels) - 1
    )
  
  widths <-
    transformed /
    sum(transformed) *
    available_width
  
  starts <- cumsum(
    c(
      left,
      head(
        widths,
        -1
      ) +
        gap
    )
  )
  
  ends <-
    starts +
    widths
  
  tibble(
    Label = labels,
    Start = starts,
    End = ends,
    Centre = (
      starts +
        ends
    ) / 2
  )
}


#Draw network panel

draw_network_panel <- function(
    site,
    time_period = NULL,
    panel_letter,
    panel_title,
    show_pollinator_labels = FALSE
) {
  
  site <- as.character(
    site
  )
  
  panel_interactions <- network_interactions %>%
    filter(
      Site == site
    )
  
  panel_flowers <- network_flowers %>%
    filter(
      Site == site
    )
  
  if (
    !is.null(
      time_period
    )
  ) {
    
    time_period <- as.character(
      time_period
    )
    
    panel_interactions <- panel_interactions %>%
      filter(
        Time.Period ==
          time_period
      )
    
    panel_flowers <- panel_flowers %>%
      filter(
        Time.Period ==
          time_period
      )
  }
  
  
  plot(
    NA,
    xlim = c(
      -0.10,
      1.30
    ),
    ylim = c(
      -0.30,
      1.34
    ),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i"
  )
  
  
  text(
    x = -0.05,
    y = 1.23,
    labels = paste0(
      panel_letter,
      ") ",
      panel_title
    ),
    adj = c(
      0,
      0.5
    ),
    cex = 1.25,
    font = 2
  )
  
  
  #Floral abundance
  
  floral_summary <- panel_flowers %>%
    group_by(
      Plant.Species
    ) %>%
    summarise(
      Floral.Abundance = mean(
        Abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    filter(
      is.finite(
        Floral.Abundance
      ),
      Floral.Abundance > 0
    )
  
  
  plants_with_visits <- unique(
    panel_interactions$Plant.Species
  )
  
  plants_present <- union(
    floral_summary$Plant.Species,
    plants_with_visits
  )
  
  plants_present <- plant_order[
    plant_order %in%
      plants_present
  ]
  
  
  plant_summary <- tibble(
    Plant.Species =
      plants_present
  ) %>%
    left_join(
      floral_summary,
      by = "Plant.Species"
    ) %>%
    mutate(
      Floral.Abundance = replace_na(
        Floral.Abundance,
        1
      )
    )
  
  
  plant_positions <- make_bar_positions(
    labels =
      plant_summary$Plant.Species,
    
    weights =
      plant_summary$Floral.Abundance,
    
    left =
      0.04,
    
    right =
      0.88,
    
    gap =
      0.015
  ) %>%
    rename(
      Plant.Species =
        Label
    )
  
  
  #Flower visitor totals
  
  pollinator_summary <- panel_interactions %>%
    group_by(
      Pollinator
    ) %>%
    summarise(
      Total.Rate = sum(
        Visit.Rate,
        na.rm = TRUE
      ),
      
      Family.Network = first(
        Family.Network
      ),
      
      .groups = "drop"
    ) %>%
    arrange(
      desc(
        Total.Rate
      ),
      Pollinator
    )
  
  
  pollinator_positions <- make_bar_positions(
    labels =
      pollinator_summary$Pollinator,
    
    weights =
      pollinator_summary$Total.Rate,
    
    left =
      0.04,
    
    right =
      0.88,
    
    gap =
      0.008
  ) %>%
    rename(
      Pollinator =
        Label
    ) %>%
    left_join(
      pollinator_summary %>%
        dplyr::select(
          Pollinator,
          Family.Network
        ),
      by = "Pollinator"
    )
  
  
  #Interaction links
  
  links <- panel_interactions %>%
    group_by(
      Plant.Species,
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
      plant_positions %>%
        dplyr::select(
          Plant.Species,
          Plant.X =
            Centre
        ),
      by = "Plant.Species"
    ) %>%
    left_join(
      pollinator_positions %>%
        dplyr::select(
          Pollinator,
          Pollinator.X =
            Centre
        ),
      by = "Pollinator"
    ) %>%
    filter(
      !is.na(
        Plant.X
      ),
      !is.na(
        Pollinator.X
      )
    )
  
  
  #Draw links
  
  if (
    nrow(
      links
    ) > 0
  ) {
    
    transformed_links <- sqrt(
      links$Weight
    )
    
    if (
      max(
        transformed_links
      ) ==
      min(
        transformed_links
      )
    ) {
      
      link_widths <- rep(
        2,
        nrow(
          links
        )
      )
      
    } else {
      
      link_widths <- 0.8 +
        3.2 *
        (
          transformed_links -
            min(
              transformed_links
            )
        ) /
        (
          max(
            transformed_links
          ) -
            min(
              transformed_links
            )
        )
    }
    
    
    for (
      i in seq_len(
        nrow(
          links
        )
      )
    ) {
      
      link_colour <- family_colours[
        links$Family.Network[i]
      ]
      
      if (
        is.na(
          link_colour
        )
      ) {
        
        link_colour <- "grey50"
      }
      
      segments(
        x0 =
          links$Plant.X[i],
        
        y0 =
          0.22,
        
        x1 =
          links$Pollinator.X[i],
        
        y1 =
          0.82,
        
        col =
          grDevices::adjustcolor(
            link_colour,
            alpha.f = 0.42
          ),
        
        lwd =
          link_widths[i],
        
        lend =
          "butt"
      )
    }
  }
  
  
  #Plant bars and labels
  
  if (
    nrow(
      plant_positions
    ) > 0
  ) {
    
    plant_colours <- ifelse(
      plant_positions$Plant.Species %in%
        restored_plants,
      
      restored_colour,
      
      other_plant_colour
    )
    
    rect(
      xleft =
        plant_positions$Start,
      
      ybottom =
        0.14,
      
      xright =
        plant_positions$End,
      
      ytop =
        0.22,
      
      col =
        plant_colours,
      
      border =
        plant_colours
    )
    
    text(
      x =
        plant_positions$Centre,
      
      y =
        0.115,
      
      labels =
        abbreviate_plant(
          plant_positions$Plant.Species
        ),
      
      srt =
        45,
      
      adj = c(
        1,
        0.5
      ),
      
      cex =
        0.95,
      
      font =
        3,
      
      xpd =
        NA
    )
  }
  
  
  #Flower visitor bars and labels
  
  if (
    nrow(
      pollinator_positions
    ) > 0
  ) {
    
    pollinator_colours <-
      family_colours[
        pollinator_positions$Family.Network
      ]
    
    pollinator_colours[
      is.na(
        pollinator_colours
      )
    ] <- "grey50"
    
    
    rect(
      xleft =
        pollinator_positions$Start,
      
      ybottom =
        0.82,
      
      xright =
        pollinator_positions$End,
      
      ytop =
        0.91,
      
      col =
        pollinator_colours,
      
      border =
        pollinator_colours
    )
    
    
    if (
      show_pollinator_labels
    ) {
      
      text(
        x =
          pollinator_positions$Centre,
        
        y =
          0.935,
        
        labels =
          abbreviate_pollinator(
            pollinator_positions$Pollinator
          ),
        
        srt =
          90,
        
        adj = c(
          0,
          0.5
        ),
        
        cex =
          0.90,
        
        col =
          pollinator_colours,
        
        xpd =
          NA
      )
    }
  }
  
  
  #Total link frequency and total link richness
  
  total_link_frequency <- sum(
    panel_interactions$Visits,
    na.rm = TRUE
  )
  
  total_link_richness <- nrow(
    links
  )
  
  
  text(
    x = 0.93,
    y = 0.66,
    
    labels = paste0(
      "Total link frequency: ",
      total_link_frequency,
      "\nTotal link richness: ",
      total_link_richness
    ),
    
    adj = c(
      0,
      1
    ),
    
    cex =
      0.82,
    
    font =
      2
  )
  
  
  #No flower visitors
  
  if (
    nrow(
      panel_interactions
    ) == 0
  ) {
    
    text(
      x = 0.47,
      y = 0.52,
      labels =
        "No flower visitors recorded",
      cex = 0.95
    )
  }
  
  
  invisible(
    NULL
  )
}


#Legend

draw_network_legend <- function(
    site
) {
  
  site_families <- network_interactions %>%
    filter(
      Site ==
        as.character(
          site
        )
    ) %>%
    distinct(
      Family.Network
    ) %>%
    arrange(
      Family.Network
    ) %>%
    pull(
      Family.Network
    )
  
  
  plot(
    NA,
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    ),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  
  
  text(
    x = 0.5,
    y = 0.88,
    labels =
      "Flower-visitor family",
    font = 2,
    cex = 1.15
  )
  
  
  if (
    length(
      site_families
    ) > 0
  ) {
    
    columns <- min(
      4,
      length(
        site_families
      )
    )
    
    rows <- ceiling(
      length(
        site_families
      ) /
        columns
    )
    
    
    x_positions <- rep(
      seq(
        0.05,
        0.78,
        length.out =
          columns
      ),
      times =
        rows
    )[
      seq_along(
        site_families
      )
    ]
    
    
    y_positions <- rep(
      seq(
        0.66,
        0.46,
        length.out =
          rows
      ),
      each =
        columns
    )[
      seq_along(
        site_families
      )
    ]
    
    
    rect(
      xleft =
        x_positions,
      
      ybottom =
        y_positions - 0.03,
      
      xright =
        x_positions + 0.025,
      
      ytop =
        y_positions + 0.03,
      
      col =
        family_colours[
          site_families
        ],
      
      border =
        family_colours[
          site_families
        ]
    )
    
    
    text(
      x =
        x_positions + 0.035,
      
      y =
        y_positions,
      
      labels =
        site_families,
      
      adj = c(
        0,
        0.5
      ),
      
      cex =
        0.90
    )
  }
  
  
  #Plant category
  
  text(
    x = 0.5,
    y = 0.27,
    labels =
      "Plant category",
    font = 2,
    cex = 1.15
  )
  
  
  rect(
    xleft = c(
      0.25,
      0.60
    ),
    
    ybottom =
      0.10,
    
    xright = c(
      0.275,
      0.625
    ),
    
    ytop =
      0.17,
    
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
    x = c(
      0.285,
      0.635
    ),
    
    y =
      0.135,
    
    labels = c(
      "Non-restored plants",
      "Restored plants"
    ),
    
    adj = c(
      0,
      0.5
    ),
    
    cex =
      0.95
  )
}


#Complete page for each site

plot_site_network <- function(
    site
) {
  
  site <- as.character(
    site
  )
  
  
  layout(
    matrix(
      1:5,
      nrow = 5,
      ncol = 1
    ),
    heights = c(
      1.45,
      1.10,
      1.10,
      1.10,
      0.65
    )
  )
  
  
  par(
    oma = c(
      0,
      0,
      2.5,
      0
    )
  )
  
  
  par(
    mar = c(
      0.6,
      0.4,
      0.8,
      0.4
    )
  )
  
  
  #Overall
  
  draw_network_panel(
    site =
      site,
    
    panel_letter =
      "a",
    
    panel_title =
      "Overall",
    
    show_pollinator_labels =
      TRUE
  )
  
  
  #Time period 1
  
  draw_network_panel(
    site =
      site,
    
    time_period =
      "1",
    
    panel_letter =
      "b",
    
    panel_title =
      "Time period 1"
  )
  
  
  #Time period 2
  
  draw_network_panel(
    site =
      site,
    
    time_period =
      "2",
    
    panel_letter =
      "c",
    
    panel_title =
      "Time period 2"
  )
  
  
  #Time period 3
  
  draw_network_panel(
    site =
      site,
    
    time_period =
      "3",
    
    panel_letter =
      "d",
    
    panel_title =
      "Time period 3"
  )
  
  
  #Legend
  
  draw_network_legend(
    site
  )
  
  
  #Site title
  
  mtext(
    paste(
      "Site",
      site
    ),
    side = 3,
    outer = TRUE,
    line = 0.5,
    cex = 1.9,
    font = 2
  )
  
  
  invisible(
    NULL
  )
}


#Check network data

table(
  network_interactions$Site,
  network_interactions$Time.Period
)


#Preview Site 1

plot_site_network(
  "1"
)


#Save Site 1 and Site 2 on separate pages

pdf(
  "Plant_flower_visitor_networks_by_site.pdf",
  width = 15,
  height = 17,
  family = "sans",
  onefile = TRUE
)

plot_site_network(
  "1"
)

plot_site_network(
  "2"
)

dev.off()


#Preview Site 1

plot_site_network(
  "1"
)


#Preview Site 2

plot_site_network(
  "2"
)
