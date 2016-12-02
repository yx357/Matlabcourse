function [diskArray] = create_diskArray(n, radius)
%-------------------------------------------------------------------------------
% usage: [diskArray] = create_diskArray(n, radius)
%
%   Create a collection of "n" random, non-overlapping, disks in a circular 
%   domain centered at (0,0).
%   
% input arguments:
% o n       number of disks.  n > 0.
% o radius      radius of the circular domain.
%
% output arguments:
% o diskArray   (n x 1) cell array of type Disk.
%
% Notes:
% o This code is part of Homework 10, CE 4121, Spring 2013.
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
    diskArray = cell(n,1);
    
    rmin = 1/sqrt(100*n)*radius;     % minimum disk radius
    rmax = 1/sqrt(    n)*radius;     % maximum disk radius
    
    for i = 1:n
        isok = false;
        
        while not( isok )
            % Create a candidate disk.
            x = (2*rand - 1)*radius;
            y = (2*rand - 1)*radius;

            vx = (rand - 0.5)*radius;
            vy = (rand - 0.5)*radius;

            r = (rmin + (rmax-rmin)*rand)*radius;
            mass = 10*pi*r^2;
            color = rand(1,3);

            candidate = Disk( [x;y], [vx;vy], r, mass, color );

            if x^2 + y^2 < (radius-r)^2
                isok = true;
                for j = 1:i-1
                    if overlap( candidate, diskArray{j} )
                        isok = false;
                        break;
                    end
                end
            end
        end
        
        diskArray{i} = candidate;
    end
end