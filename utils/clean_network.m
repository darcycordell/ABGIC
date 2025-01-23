function [L,S,T] = clean_network(L,S,T)

count = 1; delSub = [];
for i = 1:length(S)
indline1 = find(strcmp({L(:).fromSub},S(i).Name));
indline2 = find(strcmp({L(:).toSub},S(i).Name));

    if isempty(indline1) && isempty(indline2)
        
        delSub(count) = i;
        count = count+1;
    end
end

S(delSub) = [];

count = 1; delTran = [];
for i = 1:length(T)
indline1 = find(strcmp({L(:).fromSub},T(i).Sub));
indline2 = find(strcmp({L(:).toSub},T(i).Sub));

    if isempty(indline1) && isempty(indline2)
        
        delTran(count) = i;
        count = count+1;
    end
end

T(delTran) = [];