clc
clear all
close all
format long


f = 'new1';


m = dlmread(f);




 x = m(:,1);
 y = m(:,2);
 nx = m(:,3);
 ny = m(:,4);

 figure(2)
 quiver(x,y,nx,ny)
