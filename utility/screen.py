from paraview.simple import *
import sys

number = sys.argv[1]
reader = OpenDataFile("result.vti")

c_view = CreateRenderView()
c_view.OrientationAxesVisibility = 0
c_view.ViewSize = [1920, 1080]

dp = Show()
ColorBy(dp, ('POINTS', 'ImageScalars', 'Y'))
Render()

SaveScreenshot("./images/res" + str(number) + ".png" , view = c_view, magnification = 5, quality = 100)