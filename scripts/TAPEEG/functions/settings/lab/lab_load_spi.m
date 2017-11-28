function SPI = lab_load_spi

[SPI_file,SPI_filepath]=uigetfile({'*.spi;*.hdr;*.nii'},'Select spi or mri file');
if SPI_file == 0
    SPI = [];
    return
end
SPI_file = fullfile(SPI_filepath,SPI_file);

if strcmp(SPI_file(end-3:end),'.spi')
    SPI = lab_read_spi(SPI_file);
else
    SPI = lab_create_sp(SPI_file);
end
