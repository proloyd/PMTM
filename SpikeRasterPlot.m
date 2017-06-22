function [ xPoints,yPoints ] = SpikeRasterPlot( spikes )
%SpikeRaterPlot Create raster plot from binary spike data

[trials,timebins] = find(spikes);
trials = trials';
timebins = timebins';
vertSpikeHeight = 1;
SpikeStartTime = 0;
SpikePosition = 0;
halfSpikeHeight = vertSpikeHeight/2;
            
xPoints = [ timebins + SpikeStartTime;
                        timebins + SpikeStartTime;
                        NaN(size(timebins)) ];
yPoints = [ trials - halfSpikeHeight + SpikePosition;
                        trials + halfSpikeHeight + SpikePosition;
                        NaN(size(trials)) ];

 xPoints = xPoints(:);
 yPoints = yPoints(:);
 plot(xPoints,yPoints,'k');
 xlim([0 size(spikes,2)]);
 ylim([0.5 size(spikes,1)+0.5])
end

