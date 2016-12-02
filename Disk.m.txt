classdef Disk
%-------------------------------------------------------------------------------
% Disk
%
% Methods:
%   [d]    = Disk( x, v, radius, mass, color )
%   [d]    = move( d, a, dt )
%   [flag] = overlap( d1, d2 )
%            draw( d, fig )
%            dump( d, fid )
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
    properties 
        x      = [nan; nan];    % position vector of center point
        v      = [nan; nan];    % velocity vector        
        radius = nan;           % radius
        mass   = nan;           % mass
        color  = [0,0,0];       % [rgb] face color
    end
    
    methods
        %-----------------------------------------------------------------------
        % Constructor
        %-----------------------------------------------------------------------        
        function [d] = Disk( x, v, radius, mass, color )
            d.x      = x;
            d.v      = v;
            d.radius = radius;
            d.mass   = mass;
            d.color  = color;
        end
        
        %-----------------------------------------------------------------------
        % Move the disk through a time increment <dt> using the disk velocity 
        % and constant acceleration <a>.
        %-----------------------------------------------------------------------        
        function [d] = move( d, a, dt )
            d.x = d.x + d.v*dt + a*dt^2/2;
            d.v = d.v + dt*a;
        end

        %-----------------------------------------------------------------------
        % Returns true if the two disks overlap.
        %-----------------------------------------------------------------------
        function [flag] = overlap( d1, d2 )
            dx = d2.x - d1.x;
            flag = sqrt( dot(dx,dx) ) < d1.radius + d2.radius;
        end
        
        %-----------------------------------------------------------------------
        % Draw the disk on the currently active figure.
        %-----------------------------------------------------------------------        
        function draw( d )
            drawdisk( d.x(1), d.x(2), d.radius, d.color );
        end            
       
        %-----------------------------------------------------------------------
        % Dump the disk information out to the previously opened file identifed
        % by the file id <fid>.
        %-----------------------------------------------------------------------
        function dump( d, fid )
            fprintf(fid,'%10s (%5.2f %5.2f) (%5.2f %5.2f) %5.2f %5.2f [%5.3f %5.3f %5.3f] \n',...
                'DISK', d.x(1), d.x(2), d.v(1), d.v(2), d.radius, ...
                d.mass, d.color(1), d.color(2), d.color(3));
        end
    end
end