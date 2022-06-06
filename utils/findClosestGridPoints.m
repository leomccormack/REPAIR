function [idx_closest, dirs_closest, angle_diff] = findClosestGridPoints(grid_dirs_rad, target_dirs_rad)

nGrid = size(grid_dirs_rad,1);
nDirs = size(target_dirs_rad,1);
xyz_grid = unitSph2cart(grid_dirs_rad);
xyz_target = unitSph2cart(target_dirs_rad);
idx_closest = zeros(nDirs,1);
for nd=1:nDirs
    [~, idx_closest(nd)] = max(dot(xyz_grid, repmat(xyz_target(nd,:),nGrid,1), 2));
end  
dirs_closest = grid_dirs_rad(idx_closest,:);
angle_diff = acos( dot(xyz_grid(idx_closest,:), xyz_target, 2) );

end