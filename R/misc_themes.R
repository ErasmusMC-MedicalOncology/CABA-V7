# Helper theme.
theme_Job <- ggplot2::theme(
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
    axis.text.x = ggtext::element_markdown(),
    axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
    axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
    strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
    panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
    panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
    panel.border = ggplot2::element_rect(fill = NA, colour = NA),
    strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
    legend.text = ggtext::element_markdown()
)


themeTrack_Job <- theme(
    legend.position = 'right',
    legend.direction = 'horizontal',
    text = ggplot2::element_text(size = 7, family = 'Helvetica', face = 'bold'),
    axis.title.y = ggtext::element_textbox_simple(size = 6, orientation = 'left-rotated', width = NULL, halign = .5),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
    panel.grid.minor.y = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
    panel.border = ggplot2::element_rect(fill = NA, colour = NA),
    strip.background = ggplot2::element_rect(colour = 'grey20', fill = 'white'),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
)

themeAnno_Job <- theme(
    legend.position = 'right',
    axis.ticks = ggplot2::element_blank(),
    axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
    axis.text.x = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 8, family='Helvetica', face = 'bold'),
    panel.background = ggplot2::element_rect(fill = NA),
    panel.grid.major = ggplot2::element_line(NULL),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
)

colorMuts <- c(
    'Neutral' = 'white',
    'Synonymous variant' = '#A3C6BA',
    'Stop gained' = '#4897D8', # Blue
    'Missense variant' = '#000000', # Black
    'Splicing variant' ='#F9A603', # Orange
    'Multiple mutations' = 'red', # Red
    'Frameshift variant' = '#CE3B66',
    'Disruptive inframe insertion/deletion' = '#E9BAC9'
)
