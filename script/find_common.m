%�ӱ�ѡ���в��ҹ�ͬ��Լ��Ӧ

% ��ʼ���µ�Ԫ������
cell_YBplus = {};

% ��ȡ���а���'YBplus'�������Ƶ�����
indices_YBplus = find(contains(model.varNames, 'YBplus'));

% �������а���'YBplus'����
for i = indices_YBplus'
    % ���������ȫ��Ϊ1
    if all(sol_matrix(i, :) == 1)
        % ��������ӵ��µ�Ԫ����
        cell_YBplus{end+1, 1} = model.varNames{i};
    end
end


