t = [1: length(acccelY)]*0.099; #Multiply by sampling frequency
yaw_error = yaw-yaw_des;

map = imread('image.png'); %REPLACE WITH IMAGE NAME 
map2 = flipud(map);
imshow(map2);

hold on 
plot(x,y,'-')
xlabel("X Position [m]")
ylabel("Y Position [m]")
title("Path Overlayed on Image")
grid on;
axis equal;