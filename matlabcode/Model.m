classdef Model
%-------------------------------------------------------------------------------
% Model
% 
% Methods:
%   [m] = Model( radius, diskArray, acceleration )
%   [m] = proceed( m, frameRate, maxTime, fig )
%
%         Disk_Wall_Time(m, disk1)
%         Disk_Disk_Time(m, disk1, disk2)
%         Disk_Wall_Vel(m, disk1)
%         Disk_Disk_Vel(m, disk1, disk2)
%         draw( m, fig )
%         dump( m, fid )
%
% Notes:
% o This code is part of Homework 10, CE4121, Spring 2013.
% o This program was inspired by Professor Barnes, and other class peers. Thanks them.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Written by:
% Yunhan Xu
% xuxxx870@umn.edu
%
% This code was written for, and tested with MATLAB R2010b. This code 
% should function equally well on any later/earlier versions of MATLAB.
% 
% Version:
% 20130417/3:28pm
%-------------------------------------------------------------------------------
    properties
        radius       = nan;         % radius of the circular domain.
        diskArray    = cell(0,1);   % cell array of disks.
        acceleration = nan;         % global acceleration vector (2x1).
        currentTime  = nan;         % current time.
        queue        = nan;         % priority queue for events.
    end
    
    properties (Constant = true)
        EVENT_REDRAW    = 1;
        EVENT_DISK_WALL = 2;
        EVENT_DISK_DISK = 3;

        TOLERANCE       = 1e-9;
    end    
    
    methods
%-----------------------------------------------------------------------
% Constructor.
%-----------------------------------------------------------------------        
        function [m] = Model( radius, diskArray, acceleration )
            % Intialize the class properties.
              m.radius = radius;
              m.diskArray = diskArray;
              m.acceleration = acceleration;
              m.currentTime = 0;
              n_disks = length(diskArray);
            % Initialize the priority queue of events.
              m.queue = PriorityQueue(10);
            % Include the initial draw.
              ev = Event(0,1,NaN,NaN);
              m.queue = insert(m.queue, ev);
            % Include possible future disk-wall collisions.
              if n_disks >= 1
                  for i = 1:n_disks
                      dt = Disk_Wall_Time(m,i);
                      ev = Event( m.currentTime + dt, m.EVENT_DISK_WALL, i, NaN);
                      m.queue = insert(m.queue, ev);
                  end
              end
            % Include possible future disk-disk collisions.
               if n_disks >= 2
                   for i = 1:n_disks
                       for j = 1:n_disks   
                           if j ~= i
                          [dt] = Disk_Disk_Time(m, i, j);
                          ev = Event( m.currentTime + dt, m.EVENT_DISK_DISK, i, j); % ev: new event        
                          m.queue = insert(m.queue, ev);
                           end
                       end
                   end
               end     
        end

        %-----------------------------------------------------------------------        
        % Simulate through time.
        %-----------------------------------------------------------------------                
        function [m] = proceed( m, frameRate, maxTime, fig )
            n_disks = length(m.diskArray);
            
            while not( isempty(m.queue) )
                [m.queue, e] = pop(m.queue);
                % If the next event is after maxTime, put it back in the
                % queue and break out of the loop.
                  if e.time > maxTime
                      m.queue = insert(m.queue, e);
                      break
                  else    
                % Update time.
                  dt = e.time - m.currentTime;
                  m.currentTime = e.time;
                % Update the ball locations and velocities.
                  for i = 1:n_disks
                      m.diskArray{i} = move(m.diskArray{i}, m.acceleration, dt);
                  end
                  end
                % Process the event. (This is the hard part.)
                switch e.type
                    case m.EVENT_REDRAW
                        draw(m, fig);
                        dt = 1/frameRate;
                        ev = Event(m.currentTime + dt, m.EVENT_REDRAW, NaN, NaN); % ev: new event
                        m.queue = insert(m.queue, ev);
                        
                    case m.EVENT_DISK_WALL
                        %cancel events
                        m.queue = cancel(m.queue, e.id1);
                        
                        %Recalculate the velocity of the disk 
                        v = Disk_Wall_Vel(m, e.id1); 
                        m.diskArray{e.id1}.v = v;
                        %Calculate the time to wall
                        dt = Disk_Wall_Time(m,e.id1);
                        % new event
                        ev = Event(m.currentTime + dt, m.EVENT_DISK_WALL, e.id1, NaN); % ev: new event
                        m.queue = insert(m.queue, ev);
                        % Next disk_Wall hits
                        if n_disks > 1
                            for i=1:n_disks
                                if i~=e.id1
                                    [dt] = Disk_Disk_Time(m, e.id1, i);
                                    ev = Event(m.currentTime + dt, m.EVENT_DISK_DISK, e.id1, i); % ev: new event
                                    m.queue = insert(m.queue, ev);
                                end
                            end
                        end
                    
                    case m.EVENT_DISK_DISK
                        %Cancel events for involved disks, id1, and id2
                        m.queue = cancel(m.queue, e.id1);
                        m.queue = cancel(m.queue, e.id2);
                        
                        %Recalculate the velocity for each disk
                        [v1,v2] = Disk_Disk_Vel(m, e.id1, e.id2);
                        m.diskArray{e.id1}.v = v1;
                        m.diskArray{e.id2}.v = v2;
                                         
                        dt1 = Disk_Wall_Time(m, e.id1);
                        dt2 = Disk_Wall_Time(m, e.id2);
                        
                        e1= Event(m.currentTime + dt1, m.EVENT_DISK_WALL, e.id1,NaN);
                        e2= Event(m.currentTime + dt2, m.EVENT_DISK_WALL, e.id2,NaN);
                        
                        m.queue = insert(m.queue, e1);
                        m.queue = insert(m.queue, e2);
                       
                        % Next Disk_Disk hits 
                        if n_disks > 2 
                            for i=1:n_disks
                                if i~=e.id1 && i~=e.id2
                                    [dt] = Disk_Disk_Time(m, e.id1, i);
                                    ev = Event(m.currentTime + dt, m.EVENT_DISK_DISK, e.id1, i); % ev: new event
                                    m.queue = insert(m.queue, ev);
                                    
                                    [dt] = Disk_Disk_Time(m, e.id2, i);
                                    ev = Event(m.currentTime + dt, m.EVENT_DISK_DISK, e.id2, i); % ev: new event
                                    m.queue = insert(m.queue, ev);
                                end %if end
                            end %for end
                        end %if end
                end %switch end
            end %while end
        end %function end
        %%
        %-----------------------------------------------------------------------
        % Any factored routines that you deem appropriate.
        %-----------------------------------------------------------------------
        % Disk_Wall_Time function was used to find the time required to hit
        % the wall.[Provided by Dr.Barnes]
        
        function [dt] = Disk_Wall_Time(m, disk1)    
            
            d1 = m.diskArray{disk1};
            a = m.acceleration;
            
            c(1) = 0.25* dot(a,a);
            c(2) = dot(d1.v,a);
            c(3) = dot(d1.x,a) + dot(d1.v,d1.v);
            c(4) = 2*dot(d1.x,d1.v);
            c(5) = dot(d1.x,d1.x) - (m.radius-d1.radius)^2;
            t = roots(c);

            dt = inf;
               if not(isempty(t))
                  t = t(imag(t) == 0 & real(t) > m.TOLERANCE);
                  if not(isempty(t))
                     dt = min(t);
                  end
               end
               
        end    
        %------------------------------------------------------------------
        % Disk_Disk_Time function was used to find deltat for disk1 and
        % disk2 to hit together. [Provided by Dr.Barnes]
        
        function [dt] = Disk_Disk_Time(m, disk1, disk2)
            
            d1 = m.diskArray{disk1};
            d2 = m.diskArray{disk2};
            
            dx = d1.x - d2.x;
            dv = d1.v - d2.v;
            c(1) = dot(dv,dv);
            c(2) = 2*dot(dx,dv);
            c(3) = dot(dx,dx) - (d1.radius + d2.radius)^2;
            t = roots(c);
            dt = inf;
            if not( isempty(t) )
               t = t( imag(t) == 0 & real(t) > m.TOLERANCE );
               if not( isempty(t) )
                   dt = min(t);
               end
            end
        end
        %------------------------------------------------------------------
        % Disk_Wall_Vel function was used to recalculate the velocity
        % because of the wall hits.
        
        function [v] = Disk_Wall_Vel(m,disk1)
          d1 = m.diskArray{disk1};
          %create the transformation matrix
          T = 1/(d1.x(1)^2 + d1.x(2)^2)*[(d1.x(2)^2-d1.x(1)^2) -2*d1.x(1)*d1.x(2); -2*d1.x(1)*d1.x(2) (d1.x(1)^2 -d1.x(2)^2)];
          v = T*d1.v;
        end
        %------------------------------------------------------------------
        % This function was used to recalculate the velocity of each disk
        % because of the disk_disk hits [Codes provided by Dr.Barnes]
        
        function [v1,v2] = Disk_Disk_Vel(m, disk1, disk2)
             d1 = m.diskArray{disk1};
             d2 = m.diskArray{disk2};
             % Step 0
             dx = d2.x - d1.x;
             un = dx/sqrt(dx'*dx);
             ut = [-un(2); un(1)];
             % Step 1
             v1n = dot(d1.v, un);
             v1t = dot(d1.v, ut);
             v2n = dot(d2.v, un);
             v2t = dot(d2.v, ut);
             % Step 2
             v1np = ((d1.mass - d2.mass)*v1n + (2*d2.mass)*v2n) / (d1.mass + d2.mass);
             v2np = ((d2.mass - d1.mass)*v2n + (2*d1.mass)*v1n) / (d1.mass + d2.mass);
             % Step 3
             v1 = v1np*un + v1t*ut;
             v2 = v2np*un + v2t*ut;
        end
        %-----------------------------------------------------------------------
        % Compute the system kinetic energy.
        %-----------------------------------------------------------------------        
        function [pe, ke] = energy( m )
            pe = 0;
            ke = 0;
            
            for j = 1:length(m.diskArray)
                d = m.diskArray{j};
                pe = pe - d.mass * dot( m.acceleration, (d.x - m.radius*[-1;-1]) );
                ke = ke + 0.5 * d.mass * dot( d.v, d.v );
            end
        end
        
        %-----------------------------------------------------------------------
        % Draw the entire model on figure <fig>.
        %-----------------------------------------------------------------------        
        function draw( m, fig )
            % Draw the domain.
            figure(fig);
            drawdisk( 0, 0, m.radius, [1,1,1] );
            axis([-m.radius, m.radius, -m.radius, m.radius]);
            axis square;
            axis off;
            
            hold on;
            
            % Annotation.
            [pe,ke] = energy(m);
            title( sprintf('time: %.2f   energy: [%.2f, %.2f, %.2f]', m.currentTime, pe, ke, pe+ke) );            
            
            % Draw the diskArray.
            for i = 1:length(m.diskArray)
                draw( m.diskArray{i} );
            end
            
            hold off;
        end            
        
        %-----------------------------------------------------------------------
        % Dump the disk information out to the previously opened file identifed
        % by the file id <fid>.
        %-----------------------------------------------------------------------
        function dump( m, fid )
            fprintf(fid,'\n---------------------\n');
            fprintf(fid,'%10s %10.2f %5d (%10.3f,%10.3f) %10.3f %5d\n', ...
                'MODEL', m.radius, length(m.diskArray), m.acceleration(1), m.acceleration(2), ...
                m.currentTime, nactive(m.queue) );
            
            for i = 1:length(m.diskArray)
                fprintf('%10d ', i);
                dump( m.diskArray{i}, fid );
            end
            fprintf(fid,'---------------------\n');            
        end
    end
end