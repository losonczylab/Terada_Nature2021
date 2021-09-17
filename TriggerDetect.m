% Crossings = RunNoiseDetect(Signal,UpThresh,DownThresh)
%
% finds all those points where the signal crosses UpThresh in the
% upwards direction, but only for the first time after each time
% it was below DownThresh

function Crossings = TriggerDetect(Signal,UpThresh,DownThresh);

% Signal = wtas; % smoothing 5000bins
% UpThresh = 1.2*std(wtas);
% DownThresh = 1.2*median(wtas);

Signal = Signal(:);
n = length(Signal);
PreVal = [(UpThresh+DownThresh)/2; Signal(1:n-1)];

UpCrossings = find(PreVal<UpThresh & Signal>=UpThresh);
DownCrossings = find(PreVal>DownThresh & Signal<=DownThresh);

UpDownUnsort = [UpCrossings; DownCrossings];
TypeUnsort = [ones(size(UpCrossings)) ; -ones(size(DownCrossings))];

[UpDownSort Index] = sort(UpDownUnsort);
TypeSort = TypeUnsort(Index);

% we want all those times a upcrossing comes after a downcrossing
PrevType = [NaN; TypeSort(1:length(TypeSort)-1)];


Crossings.Up = UpDownSort(find(PrevType==-1 & TypeSort==1));
Crossings.Dwn = UpDownSort(find(PrevType==1 & TypeSort==-1));


if length(Crossings.Up) == 1 & length(Crossings.Dwn) == 1;
    if Crossings.Up < Crossings.Dwn;
        Crossings.Area = [ Crossings.Up Crossings.Dwn ];
    elseif Crossings.Up > Crossings.Dwn;
        Crossings.Area = [ [1 Crossings.Up]' [Crossings.Dwn length(Signal)]' ];
    end
elseif length(Crossings.Up) + length(Crossings.Dwn) == 0;
    Crossings.Area = [1 2];
elseif length(Crossings.Up) == 0 & length(Crossings.Dwn) ~= 0;
    Crossings.Area = [1 Crossings.Dwn];
elseif length(Crossings.Up) ~= 0 & length(Crossings.Dwn) == 0;
    Crossings.Area = [Crossings.Up length(Signal)];
    
elseif length(Crossings.Up) >1 | length(Crossings.Dwn) >1;
    if Crossings.Up(end) > Crossings.Dwn(end) & Crossings.Up(1) > Crossings.Dwn(1);
        Crossings.Dwn = [Crossings.Dwn;length(Signal)];
        Crossings.Up = [1; Crossings.Up];
    elseif Crossings.Up(end) > Crossings.Dwn(end) & Crossings.Up(1) < Crossings.Dwn(1);
        Crossings.Dwn = [Crossings.Dwn;length(Signal)];
    elseif Crossings.Up(end) < Crossings.Dwn(end) & Crossings.Up(1) > Crossings.Dwn(1);
        Crossings.Up = [1; Crossings.Up];
    end
    Crossings.Area = [Crossings.Up Crossings.Dwn];
end


end
