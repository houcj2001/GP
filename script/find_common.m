%从备选解中查找共同制约反应

% 初始化新的元胞数组
cell_YBplus = {};

% 获取所有包含'YBplus'的行名称的索引
indices_YBplus = find(contains(model.varNames, 'YBplus'));

% 遍历所有包含'YBplus'的行
for i = indices_YBplus'
    % 如果行数据全部为1
    if all(sol_matrix(i, :) == 1)
        % 将名称添加到新的元胞中
        cell_YBplus{end+1, 1} = model.varNames{i};
    end
end


