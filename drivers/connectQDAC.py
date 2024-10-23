import socket
import serial
import time

def send_lan(ip, port, command):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        # Connect to the device
        sock.connect((ip, port))
        
        # Send the SCPI command followed by a newline
        sock.sendall(f"{command}\n".encode())
        
        # Receive the response from the device
        response = sock.recv(1024)
        
        return response.decode()

# Example usage
ip_address = '192.168.1.253'  # Replace with your device's IP address
port = 5025                   # Replace with the correct port for SCPI commands
command = '*IDN?'              # Example SCPI command to request device identity
response = send_lan(ip_address, port, command)
print(response)



def send_usb(serial_port, baud_rate, command):
    with serial.Serial(port=serial_port, baudrate=baud_rate, timeout=1) as ser:
        # Send the SCPI command followed by a carriage return and newline
        ser.write(f"{command}\r\n".encode())
        
        # Wait for the device to respond and read the response
        response = ser.read_until().decode().strip()
        
        return response

serial_port = 'COM6'  # Replace with your device's serial port
baud_rate = 921600      # Replace with the baud rate configured for your device
command = 'sour2:volt:slew?'     # Example SCPI command to request device identit)
response = send_usb(serial_port, baud_rate, command)
print(response)

command = 'sour5:volt?'
response = send_usb(serial_port, baud_rate, command)
print(response)

command = 'sour4:volt -1.6'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -1.5'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -1.4'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -1.3'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -1.2'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -1.1'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -1'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -0.9'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -0.8'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -0.7'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt-0.6 '
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -0.5'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -0.4 '
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -0.3'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -0.2'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt -0.1'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt 0'
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt '
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt '
send_usb(serial_port, baud_rate, command)
time.sleep(100)
command = 'sour4:volt '
send_usb(serial_port, baud_rate, command)
time.sleep(100)