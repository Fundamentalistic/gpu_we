CREATE DATABASE modelling{0};
USE modelling{0};

CREATE TABLE `protein` (
   `id` BIGINT AUTO_INCREMENT,
   `res_name` VARCHAR(6) NOT NULL,
   `res_name_id` VARCHAR(6) NOT NULL,
   `res_id` BIGINT NOT NULL,
   `atom_name_id` INTEGER NOT NULL,
   `atom_id` BIGINT NOT NULL,
   `x` FLOAT, `y` FLOAT, `z` FLOAT,
   PRIMARY KEY ( `id`, `res_name_id`, `atom_name_id`, `atom_id`, `x`, `y`, `z`)
);

CREATE TABLE `water` (
   `id` BIGINT AUTO_INCREMENT,
   `res_id` BIGINT NOT NULL,
   `atom_name_id` INTEGER NOT NULL,
   `atom_id` BIGINT NOT NULL,
   `x` FLOAT, `y` FLOAT, `z` FLOAT,
   PRIMARY KEY ( `id`, `atom_name_id`, `atom_id`, `x`, `y`, `z`)
);