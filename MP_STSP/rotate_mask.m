function [final_mask] = rotate_mask(mask,robot_heading)

    final_mask = circshift(mask,-robot_heading);

end

