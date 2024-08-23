make_hex <- function() {
  imgurl <- system.file("mlmc_img.png", package = "mlmc")
  hexSticker::sticker(imgurl,
                      s_x = 1,
                      s_y = 1,
                      s_width = 0.85,
                      s_height = 0.9,
                      package="mlmc",
                      p_x = 0.6,
                      p_y = 1.45,
                      p_color = "#002147",
                      p_family = "serif",
                      p_fontface = "bold",
                      p_size = 5,
                      h_fill = "#EDEDED",
                      h_color = "#002147",
                      url = "mlmc.louisaslett.com",
                      u_size = 1.5,
                      u_family = "mono",
                      white_around_sticker = TRUE,
                      filename="inst/mlmc_hex.png",
                      dpi=600)
  usethis::use_logo("inst/mlmc_hex.png", geometry = "480x556")
}
