
# 定义 SEIR 模型函数
#' @title Ordinary Differential Equation Differentiation
#' @description 这个函数用于定义 SEIR 模型的常微分方程，
#'     描述了在传染病传播过程中易感者（S）、暴露者（E）、首次感染者（I1）、
#'     二次感染者（I2）和康复者（R）各类人群的人数变化。
#'
#' @param t 时间，自变量
#' @param y 各人群状态的初始值向量，依次为 S、E、I1、I2、R
#' @param param 参数列表，包含 mu、lamda、beta、gamma1、gamma2 和 N
#'
#' @return 返回一个列表，包含了 S、E、I1、I2 和 R 各人群状态的变化率
#' @export
#'
#' @examples
#'  \dontrun{
#' # 设置评估参数的初始值
#' times <- seq(0, 156, by = 1/7)
#' param <- c(mu = 0.000, lamda = 0.03, beta = 4, gamma1 = 0.1, gamma2 = 0.05, N = 1)
#' init <- c(S = 0.9999, E = 0.00008, I1 = 0.00001, I2 = 0, R = 0)
#'
#' # 调用常微分方程求解函数
#' result <- ODEDiff(t = times, y = init, param = param)
#' result <- as.data.frame(result)
#'
#' # 结果画图
#' seirplot <- ggplot(data = result, aes(x = time)) +
#'   geom_line(aes(y = S, color = "S"), linewidth = 1) +
#'   geom_line(aes(y = E, color = "E"), linewidth = 1) +
#'   geom_line(aes(y = I1, color = "I1"), linewidth = 1) +
#'   geom_line(aes(y = I2, color = "I2"), linewidth = 1) +
#'   geom_line(aes(y = R, color = "R"), linewidth = 1) +
#'   labs(x = "Time", y = "Ratio") +
#'   scale_color_manual(name = "SEIR",
#'                      values = c("S" = "orange", "E" = "purple", "I1" = "red", "I2" = "blue", "R" = "green")) +
#'   theme_minimal()
#'
#' # 绘制仿真结果并保存为矢量文件
#' print(seirplot)
#' ggsave(filename = "Y.ODEDiff_plot.pdf", plot = seirplot, width = 7, height = 6)
#' ggsave(filename = "Y.ODEDiff_plot.svg", plot = seirplot, width = 7, height = 6)
#' }
#'
ODEDiff <- function(t, y, param) {
  # 设定参数
  S <- y[1]
  E <- y[2]
  I1 <- y[3]
  I2 <- y[4]
  R <- y[5]
  N <- param["N"]

  # 设定参数
  beta <- param["beta"]
  mu <- param["mu"]
  gamma1 <- param["gamma1"]
  gamma2 <- param["gamma2"]
  lamda <- param["lamda"]

  # 传染病数学模型
  dSt <- mu * (N - S) - beta * S * (I1 + I2) / N
  dEt <- beta * S * (I1 + I2) / N - mu * E - lamda * E
  dI1t <- lamda * E - (mu + gamma1) * I1
  dI2t <- gamma1 * I1 - (mu + gamma2) * I2
  dRt <- gamma2 * I2 - mu * R

  # 求解结果整合成向量表达
  outcome <- c(dSt, dEt, dI1t, dI2t, dRt)

  # 返回常微分方程系统求解结果
  return(list(outcome))
}
