function cases = case_definitions()
%CASE_DEFINITIONS  The seven (C,D) parameter sets of the study, ordered by eps = D/C.
%  Amplitude A = 1, physical window t in [0, 8] (as in the original manuscript).
P = { 'Weak',                    1.0, 0.1 ;
      'Weak, stiff',             4.0, 0.5 ;
      'Moderate',                2.0, 1.0 ;
      'Moderate, stiff',        10.0, 5.0 ;
      'Strong',                  1.0, 1.0 ;
      'Strong, soft',            0.1, 0.1 ;
      'Very strong',             0.5, 2.0 };
for i = 1:size(P,1)
    cases(i).name = P{i,1}; %#ok<AGROW>
    cases(i).C = P{i,2}; cases(i).D = P{i,3};
    cases(i).eps = P{i,3}/P{i,2};
    cases(i).A = 1;
end
end
