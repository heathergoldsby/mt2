function [f a] = newfigure()
%NEWFIGURE Create a new figure, and return a handle to the figure and
%          its axis.
f = figure('color',[1 1 1]);
a = axes('Parent',f,'Position',[0.13 0.1952 0.7771 0.7229]);
box off;
hold all;

