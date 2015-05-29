#!/bin/bash

cd ../output
ffmpeg -r 12 -i graph_edges_%d.jpg    graph_edges.mp4
ffmpeg -r 12 -i depth_superpix_%d.jpg depth_superpix.mp4
ffmpeg -r 12 -i color_superpix_%d.jpg color_superpix.mp4
ffmpeg -r 12 -i color_segments_%d.jpg color_segments.mp4
ffmpeg -r 12 -i overlay_%d.jpg        overlay.mp4
ffmpeg -r 12 -i depth_in_%d.jpg       depth_in.mp4
ffmpeg -r 12 -i color_in_%d.jpg       color_in.mp4

mkdir mp4
mv *.mp4 mp4/
