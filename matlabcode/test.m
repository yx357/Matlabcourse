function test( test_case )
    close all;
   
    switch test_case
        case 1
            % One disk, with no gravity.
            fig = 1;

            acceleration = [0;0];
            radius = 1;
            frameRate = 100;
            maxTime = 10;

            diskArray = cell(1,1);
            diskArray{1} = Disk( [0;0], [1;0], [0.20], [1], [1,0,0] );

            [m] = Model( radius, diskArray, acceleration );
            [m] = proceed( m, frameRate, maxTime, fig );
    
        case 2
            % One disk, with gravity.
            fig = 1;

            acceleration = [0;-1];
            radius = 1;
            frameRate = 100;
            maxTime = 10;

            diskArray = cell(1,1);
            diskArray{1} = Disk( [0;0], [1;0], [0.20], [1], [1,0,0] );

            [m] = Model( radius, diskArray, acceleration );
            [m] = proceed( m, frameRate, maxTime, fig );

        case 3
            % Two disks, with no gravity
            fig = 1;

            acceleration = [0;0];
            radius = 1;
            frameRate = 100;
            maxTime = 10;

            diskArray = cell(2,1);
            diskArray{1} = Disk( [-0.5;0], [ 1; 0], [0.20], [1], [1,0,0] );
            diskArray{2} = Disk( [ 0.5;0], [-1; 0], [0.20], [1], [0,1,0] );
            
            [m] = Model( radius, diskArray, acceleration );
            [m] = proceed( m, frameRate, maxTime, fig );

        case 4
            % Two disks, with no gravity, and very different masses.
            fig = 1;

            acceleration = [0;0];
            radius = 1;
            frameRate = 100;
            maxTime = 10;

            diskArray = cell(2,1);
            diskArray{1} = Disk( [-0.5;0], [ 1; 0], [0.30], [1.0], [1,0,0] );
            diskArray{2} = Disk( [ 0.5;0], [-1; 0], [0.03], [0.1], [0,1,0] );
            
            [m] = Model( radius, diskArray, acceleration );
            [m] = proceed( m, frameRate, maxTime, fig );

            
        case 5
            % Two disks, with gravity.
            fig = 1;

            acceleration = [0;-1];
            radius = 1;
            frameRate = 100;
            maxTime = 10;

            diskArray = cell(2,1);
            diskArray{1} = Disk( [-0.5;0], [ 1; 0], [0.20], [1], [1,0,0] );
            diskArray{2} = Disk( [ 0.5;0], [-1; 0], [0.20], [1], [0,1,0] );
            
            [m] = Model( radius, diskArray, acceleration );
            [m] = proceed( m, frameRate, maxTime, fig );
            
        case 6
            % Three disks, with gravity.
            fig = 1;

            acceleration = [0;-1];
            radius = 1;
            frameRate = 100;
            maxTime = 10;

            diskArray = cell(3,1);
            diskArray{1} = Disk( [ 0.5;0], [-1; 0], [0.30], [30], [1,0,0] );
            diskArray{2} = Disk( [-0.5;0], [ 1; 0], [0.10], [10], [0,0,1] );
            diskArray{3} = Disk( [ 0;0.5], [ 0;-1], [0.05], [ 5], [0,1,0] );
            
            [m] = Model( radius, diskArray, acceleration );
            [m] = proceed( m, frameRate, maxTime, fig );
            
        otherwise
            % Random disks.
            fig = 1;

            n = 5 + randi(20);
            acceleration = [0;-1];
            radius = 1;
            frameRate = 100;
            maxTime = 10;

            diskArray = create_diskArray(n,radius);
            
            [m] = Model( radius, diskArray, acceleration );
            dump(m,1);

            [m] = proceed( m, frameRate, maxTime, fig );
            dump(m,1);
    end
end