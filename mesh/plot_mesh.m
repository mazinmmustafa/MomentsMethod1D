close all; clear; clc;

Data = load("bases_info.txt");
[N,~] = size(Data);

figure()
view([45 45])
axis equal
hold on
for i=1:N
  if Data(i,10)>0
    plot3([Data(i,1) Data(i,4)],
        [Data(i,2) Data(i,5)],
        [Data(i,3) Data(i,6)],
        'g','LineWidth',100)
    plot3([Data(i,4) Data(i,7)],
        [Data(i,5) Data(i,8)],
        [Data(i,6) Data(i,9)],
        'g','LineWidth',100)  
  else
    plot3([Data(i,1) Data(i,4)],
        [Data(i,2) Data(i,5)],
        [Data(i,3) Data(i,6)],
        'LineWidth',2)
    plot3([Data(i,4) Data(i,7)],
        [Data(i,5) Data(i,8)],
        [Data(i,6) Data(i,9)],
        'LineWidth',2)  
  end
  input('');     
end

hold off