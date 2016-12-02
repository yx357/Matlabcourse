function drawdisk( x, y, r, color )
%-------------------------------------------------------------------------------
% drawdisk
%   Draw a filled circle of the specified <color>, centered at (<x>,<y>) and
%   having radius <r>, on the currently active figure.
%
% Notes:
% o This code is part of Homework 10, CE4121, Spring 2013.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Written by:
%   Dr. Randal J. Barnes
%   Department of Civil Engineering
%   University of Minnesota
%
% Version:
%   7 April 2013
%-------------------------------------------------------------------------------
    persistent cost sint;
    if isempty(cost) || isempty(sint)
        theta = linspace(0, 2*pi, 100);
        cost = cos(theta);
        sint = sin(theta);
    end

    fill( x+r*cost, y+r*sint, color );
end

