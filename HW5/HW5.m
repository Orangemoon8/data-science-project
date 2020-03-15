% Homework 5
clear all; close all; clc
%% part 1
load('fashion_mnist.mat')
x_train = im2double(X_train);
x_test = im2double(X_test);
x_train = reshape(x_train,60000,28,28,1);
x_test = reshape(x_test,10000,28,28,1);
