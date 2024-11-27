# this is ADS1115 configuration codes used in the project

# Developed by Yuntao Li, OMNI Lab, Northeastern University

# For further information, please reach out to: li.yunt@northeastern.edu

import smbus
import time

# I2C Address of the ADS1115 (change if your address is different)
ADS1115_ADDRESS = 0x48

# Configurations for the ADS1115
CONFIG_REG = 0x01
MUX_DIFF_0_1 = 0x00  # Differential mode, AIN0 - AIN1
MUX_DIFF_2_3 = 0x01  # Differential mode, AIN2 - AIN3
PGA_6_144V = 0x04    # +/-6.144V range, adjust based on your input signal
DATA_RATE = 0x00     # 8 samples per second, adjust as needed

# Constants for voltage calculation (depends on PGA setting)
VOLTAGE_RANGE = 2.048  # Volts
MAX_ADC_VALUE = 32767  # 16-bit ADC

class sec_ADS1115():
    # Function to read ADC value from ADS1115 and convert to voltage
    def read_voltage(channel):
        # Create an SMBus object
        bus = smbus.SMBus(1)  # 1 for I2C bus 1 (Raspberry Pi 2 and later)

        # Set the ADS1115 configuration for the desired channel
        if channel == 0:
            config = CONFIG_REG | MUX_DIFF_0_1 | PGA_6_144V | DATA_RATE
        elif channel == 1:
            config = CONFIG_REG | MUX_DIFF_2_3 | PGA_6_144V | DATA_RATE
        else:
            raise ValueError("Invalid channel. Use 0 or 1.")

        # Write the configuration to the ADS1115
        bus.write_i2c_block_data(ADS1115_ADDRESS, 0x01, [config >> 8, config & 0xFF])

        # Read the ADC value (16-bit unsigned)
        data = bus.read_i2c_block_data(ADS1115_ADDRESS, 0x00, 2)
        adc_value = (data[0] << 8) | data[1]

        # Convert the ADC value to voltage
        voltage = adc_value * VOLTAGE_RANGE / MAX_ADC_VALUE

        return voltage
    
# try:
#     # Infinite loop to continuously read and display the input voltage values
#     while True:
#         voltage = gpt_ADS1115.read_voltage(0)
#         print(f"Channel {0}: {voltage:.2f} V")
#         time.sleep(0.1)  # Wait for 1 second before reading again

# except KeyboardInterrupt:
#     pass
