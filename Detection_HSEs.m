function [HSE_raw] = Detection_HSEs(PL_raw, mvFrame,position,tWin, Ripple_onsets, latency_threshold, binWin);

PL_nm = (PL_raw - movmean(PL_raw,mvFrame,2)) ./ std(PL_raw,0,2);
PL_nm_NoNg = mean(PL_nm,1);
PL_nm_NoNg(find(PL_nm_NoNg < 0)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FrameRest = find(smoothdata(position(2,:),'sgolay',30) <= 0.1);
maskWin = zeros(size(PL_nm_NoNg));
maskWinRest = maskWin(1,:);
maskWinRest(FrameRest) = 1;
PL_nm_rest = PL_nm_NoNg .* maskWinRest; %

% remove Frames peri-around Running Frames
myref = smoothdata(position(2,:),'sgolay',30);
FrameRun = TriggerDetect(myref,0.5,0);

if (FrameRun.Area(1,1) - tWin) < 0;
    PL_nm_rest( 1: (FrameRun.Area(1,2)+tWin) ) = 0;
elseif (FrameRun.Area(end,2) + tWin) < length(PL_nm_rest);
    PL_nm_rest( (FrameRun.Area(1,1)-tWin): (FrameRun.Area(1,2)+tWin) ) = 0;
end
for i = 2:length(FrameRun.Area(:,1))-1;
    PL_nm_rest( (FrameRun.Area(i,1)-tWin): (FrameRun.Area(i,2)+tWin) ) = 0;
end
clear i
if (FrameRun.Area(end,2) + tWin) >= length(PL_nm_rest);
    PL_nm_rest( (FrameRun.Area(end,1)-tWin): end ) = 0;
elseif (FrameRun.Area(end,2) + tWin) < length(PL_nm_rest);
    PL_nm_rest( (FrameRun.Area(end,1)-tWin): (FrameRun.Area(end,2)+tWin) ) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HSE_ref = smoothdata(PL_nm_rest,'sgolay',5);
HSE_threshold = mean(HSE_ref) + (2*std(HSE_ref));
HSE_onset = TriggerDetect(HSE_ref,HSE_threshold,0);

%%%%%%%%%%%%%
HSE_ripID = zeros([length(HSE_onset.Up),2]);
for i = 1:length(HSE_onset.Up);
    [HSE_ripID(i,1), HSE_ripID(i,2)] = min(abs(Ripple_onsets(:,1) - HSE_onset.Up(i)));
end
clear i
HSE_rip = HSE_onset.Up(find(HSE_ripID(:,1)<=latency_threshold),1);
HSE_ripID(find(HSE_ripID(:,1)<=latency_threshold),3) = 1; %

%%%%%%%%%%%%%
HSE_raw = {};
myID = find(HSE_ripID(:,3)==1);
myArea = HSE_onset.Area(find(HSE_ripID(:,3)==1),:);
for i = 1:length(myID);
    myseq_raw = PL_nm(:,myArea(i,1)-binWin:myArea(i,2)+binWin);
    HSE_raw{i} = myseq_raw;
end
clear i myseq_raw
end