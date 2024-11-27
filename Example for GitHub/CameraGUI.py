# Developed by Yuntao Li, OMNI Lab, Northeastern University

# This is a simplified demo for one camera recording and data processing.

# For further information, please reach out to: li.yunt@northeastern.edu

import tkinter as tk
from tkinter import *
from datetime import datetime
import time
import threading
import cv2
import os
import csv
import shutil
from CSI_Camera import CSI_Camera
from ADS1115 import sec_ADS1115
from PIL import Image, ImageTk  # For displaying images

### GUI ###
root = tk.Tk()
root.title("Camera and Pupillometry GUI")
root.geometry("900x500")  # Reduced size for simplicity

# Frame storage and navigation index
saved_frames = []
current_frame_index = 0

# Center Canvas for displaying images
canvas_frame = Frame(root, padx=10, pady=10, bd=2, relief=GROOVE)
canvas_frame.pack(side=TOP, expand=True, fill=BOTH)
canvas = Canvas(canvas_frame, width=320, height=240, bg="black")
canvas.pack()

# Functionality remains the same; only layout is updated
ADS1115 = sec_ADS1115()
stop_thread = False
csvfile = "mag.csv"

# Frames for organization
settings_frame = Frame(root, padx=10, pady=10, bd=2, relief=GROOVE)
settings_frame.pack(side=LEFT, fill=Y)
controls_frame = Frame(root, padx=10, pady=10, bd=2, relief=GROOVE)
controls_frame.pack(side=RIGHT, fill=Y)

# Settings Section
tk.Label(settings_frame, text="Settings", font=("Helvetica", 14, "bold")).pack()
tk.Label(settings_frame, text="Recording Time (s):").pack(anchor=W)
entry1 = tk.Entry(settings_frame, width=10)
entry1.pack(anchor=W)
entry1.insert(0, "1")

tk.Label(settings_frame, text="FPS (58 is 30fps):").pack(anchor=W)
entry2 = tk.Entry(settings_frame, width=10)
entry2.pack(anchor=W)
entry2.insert(0, "58")

tk.Label(settings_frame, text="Contrast (0.0-2.0):").pack(anchor=W)
entry3 = tk.Entry(settings_frame, width=10)
entry3.pack(anchor=W)
entry3.insert(0, "2")

tk.Label(settings_frame, text="Exposure Time (ms):").pack(anchor=W)
entry4 = tk.Entry(settings_frame, width=10)
entry4.pack(anchor=W)
entry4.insert(0, "20")

def start_recording():
# make sure only one recoridng thread is running
#if not previewing.is_set():
# start the recording task in a thread
    t1 = threading.Thread(target=run_cameras)
    t1.setDaemon(True)
    t1.start()
    global stop_thread
    stop_thread = False

# Buttons Section
tk.Label(controls_frame, text="Controls", font=("Helvetica", 14, "bold")).pack()
tk.Button(controls_frame, text="Start Preview", command=start_recording, bg="green", fg="white", width=15).pack(pady=5)
tk.Button(controls_frame, text="Stop Preview", command=lambda: stop_recording(), bg="red", fg="white", width=15).pack(pady=5)
tk.Button(controls_frame, text="Clear Data", command=lambda: clear_data(), bg="gray", fg="white", width=15).pack(pady=5)
tk.Button(controls_frame, text="Quit", command=lambda: quit_app(), bg="black", fg="white", width=15).pack(pady=5)

# Functions to display frames
def load_saved_frames(folder_path):
    global saved_frames, current_frame_index
    saved_frames = sorted([os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".tiff")])
    current_frame_index = 0
    if saved_frames:
        display_frame(current_frame_index)
    else:
        tk.messagebox.showinfo("No Frames", "No frames found in the selected folder.")

def display_frame(index):
    global saved_frames, canvas
    if 0 <= index < len(saved_frames):
        img = Image.open(saved_frames[index])
        img = img.resize((320, 240))  # Resize to fit canvas
        tk_img = ImageTk.PhotoImage(img)
        canvas.image = tk_img  # Keep a reference to prevent garbage collection
        canvas.create_image(160, 120, image=tk_img)

def show_next_frame():
    global current_frame_index
    if current_frame_index < len(saved_frames) - 1:
        current_frame_index += 1
        display_frame(current_frame_index)

def show_previous_frame():
    global current_frame_index
    if current_frame_index > 0:
        current_frame_index -= 1
        display_frame(current_frame_index)

# Add navigation buttons below canvas
nav_frame = Frame(root, padx=10, pady=10)
nav_frame.pack()

prev_button = tk.Button(nav_frame, text="Previous", command=show_previous_frame, bg="gray", fg="white", width=10)
prev_button.pack(side=LEFT, padx=5)

next_button = tk.Button(nav_frame, text="Next", command=show_next_frame, bg="gray", fg="white", width=10)
next_button.pack(side=LEFT, padx=5)

# Camera and Data Functions
def clear_data():
    with open(csvfile, "w") as file:
        file.truncate()

def quit_app():
    global stop_thread
    stop_thread = True
    root.destroy()

def stop_recording():
    global stop_thread
    stop_thread = True

# Camera defination
width = 320
height = 240

def gstreamer_pipeline(
    sensor_id=0,
    exposuretime_low = 20000000,
    exposuretime_high = 20000000,
    capture_width=width,
    capture_height=height,
    display_width=width,
    display_height=height,
    framerate=30,
    flip_method=2,
    contrast=2,
):
    return (
        "nvarguscamerasrc sensor-id=%d wbmode=7 awblock=true gainrange=\"8 8\" ispdigitalgainrange=\"4 4\" exposuretimerange=\"%d %d\" aeLock=true ! "
        "video/x-raw(memory:NVMM), format=(string)NV12, width=(int)%d, height=(int)%d, framerate=(fraction)%d/1 ! "
        "nvvidconv flip-method=%d ! "
        "videobalance contrast=%d ! "
        "video/x-raw, width=(int)%d, height=(int)%d, format=(string)NV12 ! "
        "videoconvert ! "
        "video/x-raw, format=(string)GRAY8 ! appsink"
        % (
            sensor_id,
            exposuretime_low,
            exposuretime_high,
            capture_width,
            capture_height,
            framerate,
            flip_method,
            contrast,
            display_width,
            display_height,
        )
    )

def run_cameras():
    contrast_manual = float(entry3.get())
    exposure_manual = float(entry4.get())*1000000
    window_title = "Single CSI Cameras"

    right_camera = CSI_Camera()
    right_camera.open(
        gstreamer_pipeline(
            sensor_id=0,
            exposuretime_low=exposure_manual,
            exposuretime_high=exposure_manual,
            capture_width=width,  # Decrease image resolution
            capture_height=height,  # Decrease image resolution
            flip_method=2,
            display_width=width,  # Decrease display resolution
            display_height=height,  # Decrease display resolution
            contrast=contrast_manual,
        )
    )
    right_camera.start()

    if right_camera.video_capture.isOpened():

        cv2.namedWindow(window_title, cv2.WINDOW_AUTOSIZE)

        while right_camera.video_capture.isOpened():
            _, right_image = right_camera.read()
            right_image = cv2.putText(right_image, str(datetime.now()), (50, 50), cv2.FONT_HERSHEY_SIMPLEX, 0.5,
                            (255, 0, 0), 1, cv2.LINE_AA)

            cv2.imshow(window_title, right_image)

            if stop_thread:
                break

            keyCode = cv2.waitKey(15) & 0xFF  # Adjust waitKey time for faster processing

            if keyCode == 27:
                break

            elif keyCode == ord('r'):
                print('Start recording...')
                while right_camera.video_capture.isOpened():
                    _, right_image = right_camera.read()
                    right_image = cv2.putText(right_image, str(datetime.now()), (50, 50), cv2.FONT_HERSHEY_SIMPLEX, 0.5,
                            (255, 0, 0), 1, cv2.LINE_AA)

                    cv2.imshow(window_title, right_image)

                    if stop_thread:
                        break

                    key = cv2.waitKey(15)
                    if key & 0xFF == ord('q'):
                        break

                    voltage0 = sec_ADS1115.read_voltage(0)
                    if voltage0 < 3.8:
                        print('start collecting data')
                        datestring = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
                        base_path = '/home/omni/Desktop/Cam&Motion_Monitor/%s' % (datestring)
                        right_folder = os.path.join(base_path, "right")
                        os.makedirs(right_folder, exist_ok=True)

                        T = 0
                        T_start = time.time()
                        T_total = float(entry1.get())
                        fps = float(entry2.get())
                        #fps = 58 # empirically tested that gives 30 FPS
                        interval = 1/fps
                        all_data = []

                        while T < T_total:
                            t = 0
                            t_start = time.time()
                            while t <= interval:
                                frame_start_time = time.time()
                                
                                cam_start = time.time()
                                _, right_image = right_camera.read()
                                frametime = datetime.now()
                                elapstime = round((time.time() - T_start), 2)

                                right_image = cv2.putText(right_image, str(elapstime) + " s", (50, 50),
                                                         cv2.FONT_HERSHEY_SIMPLEX, 0.8, (255, 0, 0), 2, cv2.LINE_AA)


                                # record trigger pulse
                                voltage1 = sec_ADS1115.read_voltage(1)
                                pulse = 1 if voltage1 < 3 else 0

                                cv2.imshow(window_title, right_image)

                                data = [frametime, elapstime, pulse]
                                all_data.append(data)

                                cv2.imwrite(os.path.join(right_folder, f"{frametime}.tiff"), right_image)

                                t = time.time() - t_start

                                # calculate how much time to sleep to maintain 30 fps
                                elapsed_time = time.time() - frame_start_time
                                if elapsed_time < interval:
                                    time.sleep(interval - elapsed_time)


                                if stop_thread:
                                    break

                                key = cv2.waitKey(15)
                                if key & 0xFF == ord('q'):
                                    break

                            T = time.time() - T_start
                        
                        with open(csvfile, "a") as output:
                            writer = csv.writer(output, delimiter=",", lineterminator='\n')
                            writer.writerows(all_data)
                        
                        new_path = '/home/omni/Desktop/Cam&Motion_Monitor/%s/mag.csv' %(datestring)
                        old_path = '/home/omni/Desktop/Cam&Motion_Monitor/mag.csv'
                        shutil.copy(old_path,new_path)
                        # Load frames for review
                        load_saved_frames(right_folder)

        right_camera.stop()
        right_camera.release()
        cv2.destroyAllWindows()
    else:
        print("Error: Unable to open both cameras")
        right_camera.stop()
        right_camera.release()

    

root.protocol("WM_DELETE_WINDOW", quit_app)
root.mainloop()
