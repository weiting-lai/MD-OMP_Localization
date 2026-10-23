/* MD-OMP localisation companion: one extraction pipeline, plain stereo playback. */
'use strict';
(() => {
  const $ = (id) => document.getElementById(id);
  const TAU = Math.PI * 2;
  const colours = { bg: '#F7F6F2', active: '#EFECE4', ink: '#333333', secondary: '#6B6860', muted: '#9C988E', line: '#E0DDD4', strong: '#C8C4B8' };
  const state = { scene: null, selected: 0, playback: 'source', loop: true, muted: false, volume: .85,
    top: false, cameraYaw: -.52, cameraElevation: .66 };
  const engine = { context: null, buffers: null, descriptors: [], bufferPromise: null, master: null, voices: [],
    offset: 0, startedAt: 0, playing: false, wantPlay: false, pending: null };
  let roomDirty = true, waveDirty = true, hits = [], pointer = null, lastWaveTime = -1, lastRoomTime = -1;
  const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  const qaMode = new URLSearchParams(location.search).has('qa');
  function status(message, error = false) {
    $('playbackStatus').textContent = message;
    $('playbackStatus').classList.toggle('is-error', error);
  }
  function formatTime(seconds) { return `${Math.floor(seconds/60)}:${(seconds%60).toFixed(1).padStart(4,'0')}`; }
  function coordinate(position) { return position.map((v) => (Math.abs(v)<.0005 ? 0 : v).toFixed(2)).join(', '); }
  function currentTime() {
    if (!state.scene) return 0;
    const time = engine.offset + (engine.playing ? Math.max(0,engine.context.currentTime-engine.startedAt) : 0);
    return state.loop ? time%state.scene.duration : Math.min(time,state.scene.duration);
  }
  function active(index) { return state.playback==='mix' || (state.playback==='source' && index===state.selected); }
  function activeBank() { return state.playback==='source' ? state.selected : state.playback==='mix' ? 4 : 5; }
  function selectedRecording() {
    return state.playback==='source' ? state.scene.sources[state.selected] : state.playback==='mix'
      ? { audio: state.scene.mixedAudio, waveform: state.scene.mixedWaveform }
      : { audio: state.scene.originalAudio, waveform: state.scene.originalWaveform };
  }
  function applyGates() {
    if (!engine.context) return;
    engine.voices.forEach((voice,index) => voice.gain.gain.setTargetAtTime(index===activeBank() ? 1 : 0,engine.context.currentTime,.015));
    engine.master.gain.setTargetAtTime(state.muted ? 0 : state.volume,engine.context.currentTime,.015);
  }
  async function loadAudio() {
    if (!engine.context) {
      const Context = window.AudioContext || window.webkitAudioContext;
      if (!Context) throw new Error('This browser does not support Web Audio.');
      engine.context = new Context();
      engine.master = engine.context.createGain();
      engine.master.channelCount = 2;
      engine.master.channelCountMode = 'explicit';
      engine.master.gain.value = state.muted ? 0 : state.volume;
      engine.master.connect(engine.context.destination);
    }
    await engine.context.resume();
    if (!engine.bufferPromise) {
      engine.descriptors = [...state.scene.sources.map((source) => source.audio),state.scene.mixedAudio,state.scene.originalAudio];
      status('Loading the estimated-source stereo recordings…');
      engine.bufferPromise = Promise.all(engine.descriptors.map(async (url) => {
        const response = await fetch(`${url}?v=${encodeURIComponent(state.scene.audioRevision)}`);
        if (!response.ok) throw new Error(`Cannot load ${url} (HTTP ${response.status}).`);
        const buffer = await engine.context.decodeAudioData(await response.arrayBuffer());
        if (buffer.numberOfChannels!==2 || Math.abs(buffer.duration-state.scene.duration)>1/16000) throw new Error(`Stereo alignment check failed: ${url}`);
        return buffer;
      })).then((buffers) => { engine.buffers=buffers; return buffers; }).catch((error) => { engine.bufferPromise=null; throw error; });
    }
    await engine.bufferPromise;
  }
  function destroyVoices() {
    engine.voices.forEach((voice) => { try { voice.source.stop(); } catch (_) {} voice.source.disconnect(); voice.gain.disconnect(); });
    engine.voices=[];
  }
  function createVoices() {
    destroyVoices();
    if (engine.offset>=state.scene.duration-.002) engine.offset=0;
    const when=engine.context.currentTime+.04;
    engine.voices=engine.buffers.map((buffer,index) => {
      const source=engine.context.createBufferSource(), gain=engine.context.createGain();
      source.buffer=buffer; source.loop=state.loop; source.loopEnd=state.scene.duration;
      gain.gain.value=index===activeBank() ? 1 : 0;
      source.connect(gain).connect(engine.master);
      source.start(when,engine.offset);
      return {source,gain,when,offset:engine.offset};
    });
    engine.startedAt=when; engine.playing=true; applyGates();
  }
  function start() {
    if (engine.playing) return Promise.resolve();
    engine.wantPlay=true;
    if (engine.pending) return engine.pending;
    $('playLabel').textContent='Loading';
    engine.pending=(async () => {
      try { await loadAudio(); if (engine.wantPlay) createVoices(); else status('Ready · Playback paused.'); }
      catch (error) { engine.wantPlay=false; status(`${error.message} Retry Play after checking the local server.`,true); }
      finally { engine.pending=null; renderControls(); }
    })();
    return engine.pending;
  }
  function pause() {
    engine.wantPlay=false;
    if (engine.playing) engine.offset=currentTime();
    engine.playing=false; destroyVoices();
    if (!engine.pending) status('Paused · The selected location is retained.');
    renderControls();
  }
  function selectSource(index,shouldPlay=true) {
    if (!state.scene) return;
    state.selected=index; state.playback='source'; applyGates(); renderControls();
    if (shouldPlay) void start();
  }
  function renderControls() {
    if (!state.scene) return;
    $('playButton').classList.toggle('is-playing',engine.playing);
    $('playButton').setAttribute('aria-label',engine.playing ? 'Pause audio' : 'Play audio');
    $('playLabel').textContent=engine.playing ? 'Pause' : engine.pending ? 'Loading' : 'Play';
    document.querySelectorAll('.source-row').forEach((row,index) => {
      row.classList.toggle('is-active',active(index));
      row.querySelector('.source-focus').setAttribute('aria-pressed',String(active(index)));
    });
    $('mixButton').setAttribute('aria-pressed',String(state.playback==='mix'));
    $('referenceButton').setAttribute('aria-pressed',String(state.playback==='reference'));
    $('mixState').textContent=state.playback==='mix' ? 'Listening' : 'Listen ↗';
    $('referenceState').textContent=state.playback==='reference' ? 'Listening' : 'Compare ↗';
    const label=state.playback==='source' ? `Speaker ${state.selected+1}` : state.playback==='mix' ? 'All four extracted sources' : 'Original recording · Microphone 1';
    $('activeSummary').textContent=`${label} · Stereo playback`;
    $('loopButton').classList.toggle('is-active',state.loop); $('loopButton').setAttribute('aria-pressed',String(state.loop));
    $('muteButton').classList.toggle('is-active',state.muted); $('muteButton').setAttribute('aria-pressed',String(state.muted));
    $('volumeLabel').textContent=`${Math.round(state.volume*100)}%`;
    $('timelineSlider').value=currentTime(); $('currentTime').textContent=formatTime(currentTime());
    $('downloadSource').href=selectedRecording().audio;
    $('downloadSource').textContent=`Download ${state.playback==='source' ? `Speaker ${state.selected+1}` : state.playback==='mix' ? 'all-source mix' : 'original recording'} stereo output ↓`;
    $('selectedSourceTitle').textContent=`Position · Speaker ${state.selected+1}`;
    $('selectedTrue').textContent=`(${coordinate(state.scene.sources[state.selected].truePosition)}) m`;
    $('selectedEstimate').textContent=`(${coordinate(state.scene.sources[state.selected].estimatedPosition)}) m`;
    document.querySelectorAll('#positionsBody tr').forEach((row,index) => row.classList.toggle('is-selected',active(index)));
    ['view3d','viewTop'].forEach((id,index) => { $(id).classList.toggle('is-active',index===Number(state.top)); $(id).setAttribute('aria-pressed',String(index===Number(state.top))); });
    $('sceneTag').textContent=`${state.top ? 'TOP VIEW' : '3D VIEW'} · ROOM-CENTRED COORDINATES / m`;
    if (engine.playing) status(state.muted ? 'Output muted · Playback clock continues' : `Playing · ${label} · Stereo`);
    document.documentElement.dataset.playback=state.playback;
    document.documentElement.dataset.source=String(state.selected+1);
    roomDirty=waveDirty=true; drawRoom(); drawWave();
  }

  function canvasContext(canvas) {
    const bounds = canvas.getBoundingClientRect();
    const ratio = Math.min(window.devicePixelRatio || 1, 2);
    const width = Math.round(bounds.width * ratio);
    const height = Math.round(bounds.height * ratio);
    if (canvas.width !== width || canvas.height !== height) { canvas.width = width; canvas.height = height; }
    const context = canvas.getContext('2d');
    context.setTransform(ratio, 0, 0, ratio, 0, 0);
    context.clearRect(0, 0, bounds.width, bounds.height);
    return { context, width: bounds.width, height: bounds.height };
  }

  function drawRoom() {
    if (!state.scene) return;
    const { context: ctx, width, height } = canvasContext($('roomCanvas'));
    if (!width || !height) return;
    const bounds = state.scene.room.bounds;
    const az = state.top ? 0 : state.cameraYaw;
    const elevation = state.top ? Math.PI / 2 : state.cameraElevation;
    const projectRaw = ([x, y, z]) => [x * Math.cos(az) - y * Math.sin(az),
      -(x * Math.sin(az) + y * Math.cos(az)) * Math.sin(elevation) - z * Math.cos(elevation)];
    const corners = [];
    bounds[0].forEach((x) => bounds[1].forEach((y) => bounds[2].forEach((z) => corners.push([x, y, z]))));
    const flat = corners.map(projectRaw);
    const minX = Math.min(...flat.map((p) => p[0]));
    const maxX = Math.max(...flat.map((p) => p[0]));
    const minY = Math.min(...flat.map((p) => p[1]));
    const maxY = Math.max(...flat.map((p) => p[1]));
    const scale = Math.min((width - 68) / (maxX - minX), (height - 86) / (maxY - minY));
    const project = (point) => { const [x, y] = projectRaw(point); return [width/2 + (x-(minX+maxX)/2)*scale, height/2 + (y-(minY+maxY)/2)*scale + 3]; };
    const line = (a, b, colour = colours.line, dash = [], lineWidth = 1) => {
      const p = project(a); const q = project(b);
      ctx.beginPath(); ctx.setLineDash(dash); ctx.lineWidth = lineWidth; ctx.strokeStyle = colour;
      ctx.moveTo(...p); ctx.lineTo(...q); ctx.stroke(); ctx.setLineDash([]);
    };
    const label = (text, point, dx = 0, dy = 0, colour = colours.muted, size = 10) => {
      const [x, y] = project(point); ctx.font = `${size}px "IBM Plex Mono", monospace`; ctx.fillStyle = colour; ctx.fillText(text, x+dx, y+dy);
    };
    const [x0, x1] = bounds[0]; const [y0, y1] = bounds[1]; const [floor, ceiling] = bounds[2];
    const outline = [[x0,y0,floor],[x1,y0,floor],[x1,y1,floor],[x0,y1,floor]];
    ctx.beginPath(); outline.forEach((p, i) => { const q = project(p); if (!i) ctx.moveTo(...q); else ctx.lineTo(...q); }); ctx.closePath();
    ctx.fillStyle = colours.active; ctx.globalAlpha = .37; ctx.fill(); ctx.globalAlpha = 1;
    for (let x = Math.ceil(x0); x <= x1; x += .5) line([x,y0,floor],[x,y1,floor]);
    for (let y = Math.ceil(y0); y <= y1; y += .5) line([x0,y,floor],[x1,y,floor]);
    outline.forEach((p, i) => line(p, outline[(i+1)%4], colours.strong));
    if (!state.top) {
      outline.forEach((p, i) => {
        const top = [p[0],p[1],ceiling];
        const next = outline[(i+1)%4];
        line(p,top,colours.line,[3,4]);
        line(top,[next[0],next[1],ceiling],colours.line,[3,4]);
      });
    }
    // Floor axes and scale references remain in physical metres.
    label('x / m', [x1,y0,floor], -15, 20, colours.secondary, 9);
    label('y / m', [x0,y1,floor], -30, 0, colours.secondary, 9);
    if (!state.top) label('z / m', [x0,y0,ceiling], -28, -8, colours.secondary, 9);
    label('−2', [x0,y0,floor], -16, 13, colours.muted, 8);
    label('+2.35', [x1,y0,floor], -4, 10, colours.muted, 8);

    const selected = state.scene.sources[state.selected];
    if (!state.reference) {
      state.scene.microphones.forEach((mic) => line(mic,selected.estimatedPosition,colours.strong,[2,5],.55));
    }
    state.scene.microphones.forEach((mic, index) => {
      const [x,y] = project(mic);
      ctx.strokeStyle = colours.secondary; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x-3,y); ctx.lineTo(x+3,y); ctx.moveTo(x,y-3); ctx.lineTo(x,y+3); ctx.stroke();
      if (index === 0 || index === 4 || index === 7 || index === 11) label(`M${index+1}`,mic,6,0,colours.muted,8);
    });
    hits = [];
    state.scene.sources.forEach((source, index) => {
      const measured = source.truePosition; const estimated = source.estimatedPosition;
      const [tx,ty] = project(measured); const [ex,ey] = project(estimated);
      line(measured,estimated,colours.secondary,[3,3]);
      if (!state.top) line(estimated,[estimated[0],estimated[1],floor],colours.strong,[2,4],.7);
      const on = active(index);
      const amplitudeIndex = Math.min(source.waveform.length-1, Math.floor(currentTime()/state.scene.duration * source.waveform.length));
      const amplitude = state.scene.sources[index].waveform[amplitudeIndex] || 0;
      if (on) {
        ctx.beginPath(); ctx.arc(ex,ey,14+(engine.playing && !reducedMotion ? Math.min(amplitude*35,7) : 0),0,TAU);
        ctx.strokeStyle = colours.strong; ctx.lineWidth = 1; ctx.stroke();
      }
      ctx.fillStyle = colours.bg; ctx.strokeStyle = on ? colours.ink : colours.secondary; ctx.lineWidth = 1.3;
      ctx.fillRect(tx-4,ty-4,8,8); ctx.strokeRect(tx-4,ty-4,8,8);
      ctx.font = '9px "IBM Plex Mono", monospace'; ctx.fillStyle = colours.secondary; ctx.fillText(`S${index+1}`,tx-24,ty-8);
      ctx.beginPath(); ctx.arc(ex,ey,on ? 8 : 6.5,0,TAU); ctx.fillStyle = on ? colours.ink : colours.muted; ctx.fill();
      ctx.fillStyle = colours.bg; ctx.font = '9px "IBM Plex Mono", monospace'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(String(index+1),ex,ey+.2); ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
      hits.push({ x: ex, y: ey, index });
    });
    const origin = project([0,0,0]);
    ctx.strokeStyle = colours.muted; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(origin[0]-4,origin[1]); ctx.lineTo(origin[0]+4,origin[1]); ctx.moveTo(origin[0],origin[1]-4); ctx.lineTo(origin[0],origin[1]+4); ctx.stroke();
    ctx.fillStyle = colours.muted; ctx.font = '8px IBM Plex Mono, monospace'; ctx.fillText('ORIGIN',origin[0]-14,origin[1]+17);
    roomDirty = false;
  }

  function drawWave() {
    if (!state.scene) return;
    const {context:ctx,width,height}=canvasContext($('waveformCanvas'));
    const wave=selectedRecording().waveform, count=Math.max(10,Math.floor(width/3));
    const progress=currentTime()/state.scene.duration*width, maximum=Math.max(...wave,.001);
    for (let i=0;i<count;i++) {
      const value=wave[Math.min(wave.length-1,Math.floor(i/count*wave.length))];
      const amplitude=Math.max(1,Math.min(height*.46,value/maximum*height*.46)), x=i/count*width;
      ctx.fillStyle=x<=progress ? colours.ink : colours.strong; ctx.fillRect(x,height/2-amplitude,1.3,amplitude*2);
    }
    ctx.fillStyle=colours.ink; ctx.fillRect(progress,1,1,height-2); waveDirty=false;
  }
  function animate(timestamp) {
    if (engine.playing && !state.loop && currentTime()>=state.scene.duration) { pause(); engine.offset=state.scene.duration; status('Finished · Select a position or press Play to replay.'); }
    if (engine.playing) { $('timelineSlider').value=currentTime(); $('currentTime').textContent=formatTime(currentTime()); }
    if (roomDirty || (engine.playing && timestamp-lastRoomTime>65 && !reducedMotion)) { drawRoom(); lastRoomTime=timestamp; }
    if (waveDirty || (engine.playing && timestamp-lastWaveTime>40)) { drawWave(); lastWaveTime=timestamp; }
    if (qaMode) $('demo').dataset.audioGraph=JSON.stringify(window.mdOmpDemo.snapshot());
    requestAnimationFrame(animate);
  }
  function buildSourceControls() {
    state.scene.sources.forEach((source,index) => {
      const row=document.createElement('div'); row.className='source-row'; row.dataset.source=source.id;
      const focus=document.createElement('button'); focus.type='button'; focus.className='source-focus';
      focus.setAttribute('aria-label',`Listen to Speaker ${source.id} at estimated position ${coordinate(source.estimatedPosition)} metres`);
      focus.innerHTML=`<span class="source-number">${source.id}</span><span class="source-body"><span class="source-top"><span class="source-name">Speaker ${source.id}</span><span class="source-error">${(source.errorMetres*100).toFixed(1)} cm error</span></span><span class="source-coordinates">(${coordinate(source.estimatedPosition)}) m</span></span>`;
      focus.addEventListener('click',() => selectSource(index)); row.append(focus); $('sourceList').append(row);
      const tr=document.createElement('tr');
      [source.label,`(${coordinate(source.truePosition)})`,`(${coordinate(source.estimatedPosition)})`,(source.errorMetres*100).toFixed(2)].forEach((text,i) => { const td=document.createElement('td'); td.textContent=text; if(i===3)td.className='num-col'; tr.append(td); });
      const td=document.createElement('td'), listen=document.createElement('button');
      listen.type='button'; listen.textContent=`Listen S${source.id}`; listen.setAttribute('aria-label',`Listen to estimated Speaker ${source.id}`);
      listen.addEventListener('click',() => { selectSource(index); $('demo').scrollIntoView({behavior:reducedMotion ? 'instant' : 'smooth',block:'start'}); });
      td.append(listen); tr.append(td); $('positionsBody').append(tr);
    });
  }
  function bindEvents() {
    $('playButton').addEventListener('click',() => { if(engine.playing || engine.wantPlay)pause(); else void start(); });
    $('loopButton').addEventListener('click',() => { const offset=currentTime(); state.loop=!state.loop; engine.offset=offset; if(engine.playing)createVoices(); renderControls(); });
    $('muteButton').addEventListener('click',() => { state.muted=!state.muted; applyGates(); renderControls(); });
    $('volumeSlider').addEventListener('input',(event) => { state.volume=Number(event.target.value); applyGates(); $('volumeLabel').textContent=`${Math.round(state.volume*100)}%`; });
    $('timelineSlider').addEventListener('input',(event) => {
      const playing=engine.playing || engine.wantPlay, time=Number(event.target.value);
      pause(); engine.offset=time; renderControls(); if(playing)void start();
    });
    $('mixButton').addEventListener('click',() => { state.playback='mix'; applyGates(); renderControls(); void start(); });
    $('referenceButton').addEventListener('click',() => { state.playback=state.playback==='reference' ? 'source' : 'reference'; applyGates(); renderControls(); void start(); });

    $('view3d').addEventListener('click', () => { state.top = false; renderControls(); });
    $('viewTop').addEventListener('click', () => { state.top = true; renderControls(); });
    $('resetView').addEventListener('click', () => { state.top = false; state.cameraYaw = -.52; state.cameraElevation = .66; renderControls(); });
    const canvas = $('roomCanvas');
    canvas.addEventListener('pointerdown', (event) => {
      pointer = { x:event.clientX,y:event.clientY,lastX:event.clientX,lastY:event.clientY,moved:false,id:event.pointerId };
      canvas.setPointerCapture(event.pointerId);
    });
    canvas.addEventListener('pointermove', (event) => {
      const rect = canvas.getBoundingClientRect();
      if (!pointer) {
        canvas.style.cursor = hits.some((hit) => Math.hypot(event.clientX-rect.left-hit.x,event.clientY-rect.top-hit.y)<17) ? 'pointer' : 'grab';
        return;
      }
      if (Math.hypot(event.clientX-pointer.x,event.clientY-pointer.y)>4) pointer.moved = true;
      if (pointer.moved) {
        if (state.top) { state.top = false; renderControls(); }
        state.cameraYaw += (event.clientX-pointer.lastX)*.008;
        state.cameraElevation = Math.max(.2,Math.min(1.3,state.cameraElevation+(event.clientY-pointer.lastY)*.006));
        roomDirty = true;
      }
      pointer.lastX = event.clientX; pointer.lastY = event.clientY;
    });
    canvas.addEventListener('pointerup', (event) => {
      if (!pointer) return;
      if (!pointer.moved) {
        const rect = canvas.getBoundingClientRect();
        const nearest = [...hits].sort((a,b) => Math.hypot(event.clientX-rect.left-a.x,event.clientY-rect.top-a.y)-Math.hypot(event.clientX-rect.left-b.x,event.clientY-rect.top-b.y))[0];
        if (nearest && Math.hypot(event.clientX-rect.left-nearest.x,event.clientY-rect.top-nearest.y)<19) selectSource(nearest.index);
      }
      pointer = null;
    });
    canvas.addEventListener('pointercancel', () => { pointer = null; });
    canvas.addEventListener('keydown', (event) => {
      if (event.key === 'ArrowLeft' || event.key === 'ArrowRight') { event.preventDefault(); state.cameraYaw += event.key === 'ArrowLeft' ? -.1 : .1; state.top = false; renderControls(); }
      if (event.key === 'ArrowUp' || event.key === 'ArrowDown') { event.preventDefault(); state.cameraElevation = Math.max(.2,Math.min(1.3,state.cameraElevation+(event.key === 'ArrowUp' ? .08 : -.08))); roomDirty = true; }
    });
    document.addEventListener('keydown', (event) => {
      if (event.ctrlKey || event.altKey || event.metaKey || event.target.closest('input, textarea, select, [contenteditable="true"]')) return;
      if (/^[1-4]$/.test(event.key)) { event.preventDefault(); selectSource(Number(event.key)-1); }
      if (event.code === 'Space' && !event.target.closest('button,a')) { event.preventDefault(); if (engine.playing || engine.wantPlay) pause(); else void start(); }
    });

    $('copyCitation').addEventListener('click',async () => {
      const citation='@article{mdomp,\n  author = {Lai, Wei-Ting and Birnie, Lachlan and Abhayapala, Thushara and Bastine, Amy and Samarasinghe, Prasanga},\n  title = {Sound Source Localization using Multi-Dictionary Orthogonal Matching Pursuit in Reverberant Environments},\n  journal = {IEEE Transactions on Audio, Speech and Language Processing},\n  year = {2025},\n  doi = {10.1109/TASLPRO.2025.3624966}\n}';
      try { await navigator.clipboard.writeText(citation); $('copyCitation').textContent='Copied'; setTimeout(() => {$('copyCitation').textContent='Copy BibTeX';},2200); }
      catch (_) { $('copyCitation').textContent='Citation in README'; status('The full citation is available in the README.'); }
    });
    new ResizeObserver(() => {roomDirty=waveDirty=true;}).observe($('demo'));
    document.fonts.ready.then(() => {roomDirty=waveDirty=true;});
    window.addEventListener('pagehide',() => {pause(); if(engine.context)void engine.context.suspend();});
  }
  // Read-only diagnostics for checking the actual shared-clock stereo graph.
  window.mdOmpDemo=Object.freeze({snapshot:() => ({playback:state.playback,selectedSource:state.selected+1,
    time:currentTime(),playing:engine.playing,muted:state.muted,volume:state.volume,
    buffers:engine.buffers?.map((b) => ({channels:b.numberOfChannels,duration:b.duration})) || [],
    urls:engine.descriptors.slice(),activeBank:activeBank(),masterGain:engine.master?.gain.value,
    voices:engine.voices.map((v) => ({when:v.when,offset:v.offset,gain:v.gain.gain.value})),
    markerHits:hits.map((h) => ({...h})),camera:{top:state.top,yaw:state.cameraYaw,elevation:state.cameraElevation}})});
  async function initialise() {
    try {
      const response=await fetch('data/scene.json',{cache:'no-cache'});
      if(!response.ok)throw new Error(`Cannot load localisation results (HTTP ${response.status}).`);
      const scene=await response.json();
      if(scene.schemaVersion!==2 || scene.sources?.length!==4 || scene.microphones?.length!==14 || scene.channels!==2)throw new Error('Invalid localisation scene.');
      state.scene=scene;
      $('meanError').textContent=$('tableMeanError').textContent=`${(scene.meanErrorMetres*100).toFixed(1)} cm`;
      $('audioDetails').textContent=`Stereo · ${scene.sampleRate/1000} kHz · ${scene.duration.toFixed(2)} s`;
      $('duration').textContent=formatTime(scene.duration); $('timelineSlider').max=scene.duration;
      $('listeningStoi').textContent=scene.listening.meanStoi.toFixed(3);
      buildSourceControls(); bindEvents(); renderControls(); $('playButton').disabled=false;
      status('Ready · Select an MD-OMP estimate to listen.'); requestAnimationFrame(animate);
    } catch(error) { $('playLabel').textContent='Unavailable'; status(`${error.message} Serve this folder over HTTP; see README.`,true); }
  }
  void initialise();
})();
