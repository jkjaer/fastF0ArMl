function plotObjective(Obj, pitchOrder, arOrder)
    if isempty(Obj.objValNoPitch)
        error('Please run the estimator (Obj.estimate) first!'); 
    end
    if nargin == 1
        pitchOrder = Obj.estPitchOrder;
        arOrder = Obj.estArOrder;
    end
    if pitchOrder == 0
        disp('The pitch order is zero so nothing to plot :(');
    else
        idx = (Obj.dftRange(pitchOrder,1):Obj.dftRange(pitchOrder,2))+1;
        pitchGrid = Obj.fullPitchGrid(idx);
        plot(pitchGrid, Obj.objValPitch(arOrder+1,1:length(idx), ...
            pitchOrder), 'lineWidth', 2);
        xlabel('Pitch [Hz]');
        ylabel('Objective value')
    end
end