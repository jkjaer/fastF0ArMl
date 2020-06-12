function reportProblem(Obj)
    text = ...
        'Please describe the problem as accurately as possible: ';
    errorDescription = inputdlg(text, 'Report Error', [10, 100]);
    matFileName = ['FastF0ArMlProblem_', ...
        datestr(now,'yyyymmdd_HHMMSS'), '.mat'];
    FastF0ArMl = Obj;
    save(matFileName, 'errorDescription', 'FastF0ArMl');
    disp('A MAT-file is now stored in the current directory.');
    disp('Please send it to jkn@create.aau.dk.');
    disp('Thanks!');
end
