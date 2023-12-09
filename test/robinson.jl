using Plots

latitudes = 0.0:5.0:90.0
X = [
    1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600, 0.9427, 0.9216, 0.8962,
    0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322
]
Y = [
    0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958, 0.5571, 
    0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000
]

canvasX = plot(latitudes, X, markershape=:diamond, label="X")
canvasY = plot(latitudes, Y, markershape=:diamond, label="Y")

lats = 0.0:0.01:90.0

splineX = CubicSpline(latitudes, X)
splineY = CubicSpline(latitudes, Y)
X_cubic_splines = [splineX(x) for x in lats]
Y_cubic_splines = [splineY(x) for x in lats]

linX = LinearInterpolater(latitudes, X)
linY = LinearInterpolater(latitudes, Y)
X_lin = [linX(x) for x in lats]
Y_lin = [linY(x) for x in lats]

X_lagrange = [lagrange_polynomial(latitudes, X, x) for x in lats]
Y_lagrange = [lagrange_polynomial(latitudes, Y, x) for x in lats]

plot!(canvasX, lats, X_cubic_splines, label="cubic spline")
plot!(canvasX, lats, X_lin, label="linear")
plot!(canvasX, lats, X_lagrange, label="Lagrange")

plot!(canvasY, lats, Y_cubic_splines, label="cubic spline")
plot!(canvasY, lats, Y_lin, label="linear")
plot!(canvasY, lats, Y_lagrange, label="Lagrange")

plot(canvasX, canvasY)

#canvas_diffX = plot(lats, X_cubic_splines - X_lin)
#canvas_diffY = plot(lats, Y_cubic_splines - Y_lin)
#plot(canvas_diffX, canvas_diffY)