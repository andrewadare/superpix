import cv2

dep = cv2.imread("/Users/adare/kinect-data/clutter/depth_57.png")
rgb = cv2.imread("/Users/adare/kinect-data/clutter/color_57.jpg")

# rgbgray = cv2.cvtColor(rgb, cv2.COLOR_BGR2GRAY)
# # depgray = cv2.cvtColor(dep, cv2.COLOR_BGR2GRAY)
# blend = cv2.addWeighted(rgbgray, 0.7, dep, 0.3, 0)
# cv2.rectangle(img, pt1, pt2, color[, thickness[, lineType[, shift]]]) → None

tl = (100, 90)
br = (800, 380)
cv2.rectangle(rgb, tl, br, (255,0,0), 2)
cv2.rectangle(dep, tl, br, (0,255,0), 2)

# rgb = cv2.rectangle(rgb,(384,0),(320,128),(0,255,0),3)

# import numpy as np
# img = np.zeros((512,512,3), np.uint8)
# # Draw a diagonal blue line with thickness of 5 px
# cv2.line(img,(0,0),(511,511),(255,0,0),5)
# cv2.rectangle(img,(384,0),(510,128),(0,255,0),3)

cv2.namedWindow("foo")
cv2.imshow("foo", dep)
cv2.waitKey(0)
cv2.destroyAllWindows()
