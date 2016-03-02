function out = struct_merge( s1, s2 )
%STRUCT_MERGE Merges s1 and s2 structs, if a field exists in both, s2
%overwrites s1's value
%   Detailed explanation goes here


  S1MapDefault = containers.Map(fieldnames(s1),struct2cell(s1)); %Create unique keys
  S2MapInput = containers.Map(fieldnames(s2),struct2cell(s2)); 
  SMap = [S1MapDefault; S2MapInput];
  out = cell2struct(values(SMap),keys(SMap),2);

end

